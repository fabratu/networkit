#pragma once

#include <d2hb/block.h>
#include <d2hb/graph.h>
#include <d2hb/hash_tools.h>

#include <omp.h>

#include <algorithm>
#include <assert.h>
#include <atomic>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace d2hb {

template <typename EdgeIt>
std::tuple<EdgeIt, EdgeIt> thread_batch(EdgeIt batch_begin, EdgeIt batch_end,
                                        unsigned int thread_count, unsigned int thread_id) {
    size_t const batch_size = std::distance(batch_begin, batch_end);
    size_t const elements = batch_size / thread_count;
    if (elements == 0 || thread_count >= batch_size) {
        bool const first_thread = thread_id == 0;
        if (first_thread) {
            return {batch_begin, batch_end};
        } else {
            return {batch_end, batch_end};
        }
    }

    size_t const position = thread_id * elements;

    EdgeIt start = std::min(batch_begin + position, batch_end);

    bool const last_thread = thread_id == (thread_count - 1);
    EdgeIt end = last_thread ? batch_end : std::min(start + elements, batch_end);

    if (start != end) {
        Vertex const predecessor =
            (start == batch_begin) ? invalidVertex() : std::prev(start, 1)->source;

        while (start != end && predecessor == start->source) {
            std::advance(start, 1);
        }

        if (start != end) {
            for (Vertex successor = (end == batch_end) ? invalidVertex() : end->source;
                 end != batch_end && successor == (end - 1)->source && end->source != predecessor;
                 successor = end->source) {
                std::advance(end, 1);
            }
        }
    }

    return {start, end};
}

// Use this thread batch generator if you use the batch together with an apply
// function that does not insert new edges and therefore the exclusivity of
// vertex u to thread x must not be guaranteed.
template <typename EdgeIt>
std::tuple<EdgeIt, EdgeIt> non_exclusive_thread_batch(EdgeIt batch_begin, EdgeIt batch_end,
                                                      unsigned int thread_count,
                                                      unsigned int thread_id) {
    size_t const batch_size = std::distance(batch_begin, batch_end);
    size_t const elements = batch_size / thread_count;
    if (elements == 0 || thread_count >= batch_size) {
        bool const first_thread = thread_id == 0;
        if (first_thread) {
            return {batch_begin, batch_end};
        } else {
            return {batch_end, batch_end};
        }
    }

    bool const first_thread = thread_id == 0;
    bool const last_thread = thread_id == (thread_count - 1);

    if (elements == 0 || elements == batch_size) {
        if (first_thread) {
            return {batch_begin, batch_end};
        } else {
            return {batch_end, batch_end};
        }
    }

    EdgeIt start = batch_begin;
    if (!first_thread) {
        size_t const offset = thread_id * elements;
        std::advance(start, offset);
    }

    EdgeIt end = start;
    if (last_thread) {
        end = batch_end;
    } else {
        std::advance(end, elements);
    }

    return {start, end};
}

// TODO: The following fast modulo reduction can not be used since we have 64
// bit keys. Is there another fast modulo reduction for 64bit?
//
// Since the hashed key has high entropy we can't deduce any pattern thus we
// leave the modulo reduction to the compiler.
inline uint64_t key_to_thread(uint64_t k, uint32_t t_count) { return hash64(k) % t_count; };

template <typename T, typename ConstIt> class HypersparseBatchParallelizer {
  public:
    template <typename F> void apply_non_exclusive(ConstIt begin, ConstIt end, F func) const {
        int const t_count = omp_get_max_threads();

#pragma omp parallel shared(begin, end)
        {
            std::tuple<ConstIt, ConstIt> const local_batch =
                non_exclusive_thread_batch(begin, end, t_count, omp_get_thread_num());
            ConstIt local_begin = ConstIt(std::get<0>(local_batch));
            ConstIt local_end = ConstIt(std::get<1>(local_batch));
            func(local_begin, local_end);
        }
    }

    template <typename F> void apply_exclusive(ConstIt begin, ConstIt end, F func) const {
        int const t_count = omp_get_max_threads();

#pragma omp parallel shared(begin, end)
        {
            std::tuple<ConstIt, ConstIt> const local_batch =
                thread_batch(begin, end, t_count, omp_get_thread_num());
            ConstIt local_begin = ConstIt(std::get<0>(local_batch));
            ConstIt local_end = ConstIt(std::get<1>(local_batch));
            func(local_begin, local_end);
        }
    }

    template <typename Iterator, typename K, typename Cmp, typename F>
    void sort_and_apply_exclusive(Iterator begin, Iterator end, K key, Cmp cmp, F func) const {
        int const t_count = omp_get_max_threads();

        std::sort(begin, end, cmp);

#pragma omp parallel shared(begin, end)
        {
            std::tuple<Iterator, Iterator> const local_batch =
                thread_batch(begin, end, t_count, omp_get_thread_num());
            ConstIt local_begin = ConstIt(std::get<0>(local_batch));
            ConstIt local_end = ConstIt(std::get<1>(local_batch));
            func(local_begin, local_end);
        }
    }
};

template <typename T> class BatchParallelizer {
  public:
    template <typename Iterator, typename K, typename Cmp, typename F>
    void operator()(Iterator begin, Iterator end, K key, Cmp cmp, F func) {
        int const t_count = omp_get_max_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            for (auto it = begin; it != end; ++it)
                func(*it);
            return;
        }

        std::sort(begin, end, cmp);
#pragma omp parallel shared(begin, end)
        {
            std::tuple<Iterator, Iterator> local_batch =
                thread_batch(begin, end, t_count, omp_get_thread_num());
            for (auto it = std::get<0>(local_batch); it != std::get<1>(local_batch); ++it) {
                func(*it);
            }
        }
    }

    template <typename Iterator, typename K, typename Cmp, typename F>
    void apply(Iterator begin, Iterator end, K key, Cmp cmp, F func) {
        int const t_count = omp_get_max_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            for (auto it = begin; it != end; ++it)
                func(*it);
            return;
        }

        m_batch_counts.resize((t_count + 1) * t_count);
        m_batch_slots.resize(n);

#pragma omp parallel num_threads(t_count)
        {
            auto t = omp_get_thread_num();
            assert(omp_get_num_threads() == t_count);

            auto counts_of_thread = [&](int ct) -> Vertex* {
                return &m_batch_counts[ct * (t_count + 1)];
            };

            auto n_per_thread = n / t_count;
            auto i_begin = t * n_per_thread;
            auto i_end = i_begin + n_per_thread;
            if (t == t_count - 1)
                i_end = n;

            // First, perform a local counting sort to sort updates according to associated threads.
            auto t_counts = counts_of_thread(t);
            for (int at = 0; at < t_count; ++at)
                t_counts[at] = 0;

            for (size_t i = i_begin; i < i_end; ++i) {
                auto k = key(*(begin + i));
                auto at = key_to_thread(k, t_count);
                ++t_counts[at];
            }

            uint64_t psum = 0;
            for (int at = 0; at < t_count; ++at) {
                psum += t_counts[at];
                t_counts[at] = i_begin + psum;
            }
            assert(i_begin + psum == i_end);
            t_counts[t_count] = i_end;

            for (size_t irev = i_end; irev > i_begin; --irev) {
                // Iterating in reverse ensures that the sort is stable;
                // this yields a better memory access pattern when performing random access later.
                auto i = irev - 1;
                auto k = key(*(begin + i));
                auto at = key_to_thread(k, t_count);
                m_batch_slots[--t_counts[at]] = i;
            }

            // Now, let each thread collect its updates.

#pragma omp barrier

            uint64_t local_count = 0;
            for (int ot = 0; ot < t_count; ++ot) {
                auto ot_counts = counts_of_thread(ot);
                auto j_begin = ot_counts[t];
                auto j_end = ot_counts[t + 1];
                for (size_t j = j_begin; j < j_end; ++j) {
                    auto i = m_batch_slots[j];
                    auto edge = *(begin + i);
                    func(*(begin + i));
                }
                local_count += j_end - j_begin;
            }
        }
    }

    template <typename Iterator, typename K> void distribute(Iterator begin, Iterator end, K key) {
        int const t_count = omp_get_num_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
            return;
        }

#pragma omp single
        {
            m_batch_counts.resize((t_count + 1) * t_count);
            m_batch_slots.resize(n);
        }

        auto t = omp_get_thread_num();

        auto counts_of_thread = [&](int ct) -> uint64_t* {
            return &m_batch_counts[ct * (t_count + 1)];
        };

        auto n_per_thread = n / t_count;
        auto i_begin = t * n_per_thread;
        auto i_end = i_begin + n_per_thread;
        if (t == t_count - 1)
            i_end = n;

        // First, perform a local counting sort to sort updates according to associated threads.

        auto t_counts = counts_of_thread(t);
        for (int at = 0; at < t_count; ++at)
            t_counts[at] = 0;

        for (size_t i = i_begin; i < i_end; ++i) {
            auto k = key(*(begin + i));
            auto at = key_to_thread(k, t_count);
            ++t_counts[at];
        }

        uint64_t psum = 0;
        for (int at = 0; at < t_count; ++at) {
            psum += t_counts[at];
            t_counts[at] = i_begin + psum;
        }
        assert(i_begin + psum == i_end);
        t_counts[t_count] = i_end;

        for (size_t irev = i_end; irev > i_begin; --irev) {
            // Iterating in reverse ensures that the sort is stable;
            // this yields a better memory access pattern when performing random access later.
            auto i = irev - 1;
            auto k = key(*(begin + i));
            auto at = key_to_thread(k, t_count);
            m_batch_slots[--t_counts[at]] = i;
        }

        // Now, let each thread collect its updates.
#pragma omp barrier
    }

    template <typename Iterator, typename F> void map(Iterator begin, Iterator end, F func) {
        int const t_count = omp_get_num_threads();
        size_t const n = end - begin;
        if (t_count == 1 || n < t_count) {
#pragma omp master
            {
                for (auto it = begin; it != end; ++it) {
                    T elem = *it;
                    func(elem);
                }
            }
            return;
        }

        auto t = omp_get_thread_num();

        auto counts_of_thread = [&](int ct) -> uint64_t* {
            return &m_batch_counts[ct * (t_count + 1)];
        };

        uint64_t local_count = 0;
        for (int ot = 0; ot < t_count; ++ot) {
            auto ot_counts = counts_of_thread(ot);
            auto j_begin = ot_counts[t];
            auto j_end = ot_counts[t + 1];
            for (size_t j = j_begin; j < j_end; ++j) {
                auto i = m_batch_slots[j];
                auto edge = *(begin + i);
                func(*(begin + i));
            }
            local_count += j_end - j_begin;
        }
    }

  private:
    std::vector<uint64_t> m_batch_counts;
    std::vector<uint64_t> m_batch_slots;
};

} // namespace d2hb
