#ifndef NETWORKIT_AUXILIARY_PARALLEL_HASH_MAP_HPP_
#define NETWORKIT_AUXILIARY_PARALLEL_HASH_MAP_HPP_

#include <omp.h>

#include <limits>
#include <memory>
#include <random>
#include <vector>

#include <networkit/auxiliary/AtomicUtils.hpp>
#include <networkit/auxiliary/HashUtils.hpp>

namespace Aux {

// An object of this class will handle a dynamic hashtable and offers factory
// methods to acquire hashtable resources.
class ParallelHashMap {
public:
    friend void swap(ParallelHashMap &p, ParallelHashMap &q) noexcept {
        using std::swap;

        swap(p.m_source, q.m_source);
        swap(p.m_target, q.m_target);
        swapAtomicsNonAtomically<uint32_t>(p.m_busy_bitset, q.m_busy_bitset);
        swapAtomicsNonAtomically<bool>(p.m_request_growth, q.m_request_growth);
    }

    class HTHandle;
    class HTAtomic128;
    class HTSyncData;

    static constexpr uint64_t ht_invalid_key = std::numeric_limits<uint64_t>::max();
    static constexpr uint64_t ht_invalid_value = std::numeric_limits<uint64_t>::max() - 1;

    // This value shall represent some storage capacity that can be handled within a
    // CPU cache (e.g. private MLC).
    static constexpr size_t const ht_begin_capacity = 4096u;

    static constexpr uint64_t ht_key_space = sizeof(uint64_t) * 8;

    ParallelHashMap();

    // begin_capacity: Must be multiple of 2.
    ParallelHashMap(size_t const begin_capacity);

    ~ParallelHashMap() = default;

    ParallelHashMap(ParallelHashMap const &other);

    ParallelHashMap(ParallelHashMap &&other);

    ParallelHashMap &operator=(ParallelHashMap other) {
        swap(*this, other);
        return *this;
    }

    std::unique_ptr<HTHandle> makeHandle();

    HTAtomic128 const *const currentTable() const;

    uint32_t randomThreadRange();

private:
    std::unique_ptr<HTAtomic128> m_source;
    std::unique_ptr<HTAtomic128> m_target;
    std::atomic_uint32_t m_busy_bitset{0u};
    std::atomic_bool m_request_growth{false};

    // Used to determine how often to check for the hashtables fill factor for
    // every 1 to p inserts. We're keeping this an internal random number
    // generator engine. However, this could also be injected when constructing.
    std::mt19937 m_generator;
};

// This is the Hashtable Atomic 128 (HTAtomic128) class which stores (key,
// value) pairs, both of 64 bit and operates on them atomically using DWCAS
// operations.
//
// The basis of this implementation is the "Concurrent Hash Tables: Fast and
// general(?)!" from Maier et al. (2016).
//
// Some interfaces are not needed for our case such as insertOrUpdate() or
// update(). We expect the key's never to be updated during runtime: the handle
// for vertex u will never change during the execution of the program.
class ParallelHashMap::HTAtomic128 {
public:
    // In order to use DWCAS we must align our Cell within a 16 bytes window.
    struct alignas(16) Cell {
        uint64_t key{ParallelHashMap::ht_invalid_key};
        uint64_t value{ParallelHashMap::ht_invalid_value};
    };

    using Cells = std::vector<Cell>;

    HTAtomic128();
    HTAtomic128(size_t const capacity);
    ~HTAtomic128() = default;
    HTAtomic128(HTAtomic128 const &other);
    HTAtomic128(HTAtomic128 &&other);

    HTAtomic128 &operator=(HTAtomic128 other);

    // Thanks to GManNickG on Stackoverflow explaining the copy-swap idiom:
    // https://stackoverflow.com/a/3279550
    friend void swap(HTAtomic128 &a, HTAtomic128 &b) {
        using std::swap;

        swap(a.m_cells, b.m_cells);
        swap(a.m_capacity, b.m_capacity);
        swap(a.m_log_capacity, b.m_log_capacity);
        swapAtomicsNonAtomically(a.m_global_occupancy, b.m_global_occupancy);
    }

    // Returns false if key already exists, and true if (key, value) was
    // successfully inserted.
    bool insert(uint64_t const key, uint64_t const value);

    // Returns the actual value of the cell if key is present, ht_invalid_key
    // otherwise.
    uint64_t find(uint64_t const key) const;

    // Returns false if key does not exist, true if (key, value) was
    // successfully updated.
    bool update(uint64_t const key, uint64_t const value);

    struct Iterator {
        using iterator_category = std::input_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Cell const;
        using pointer = Cell const *;
        using reference = Cell const &;

        Iterator(Cells const &cells);

        friend void swap(Iterator &a, Iterator &b) {
            using std::swap;

            swap(a.m_cells_it, b.m_cells_it);
            swap(a.m_cells_end, b.m_cells_end);
        }

        reference operator*() { return *m_cells_it; }
        pointer operator->();

        // Prefix increment operator
        Iterator &operator++();

        // Postfix increment operator
        Iterator operator++(int);

        bool operator==(const HTAtomic128::Iterator &b) const {
            return this->m_cells_it == b.m_cells_it;
        }

        bool operator!=(const HTAtomic128::Iterator &b) const { return !(*this == b); }

        void invalidate();

    private:
        void next();

    private:
        Cells::const_iterator m_cells_it;
        Cells::const_iterator m_cells_end;
    };

    Iterator begin() const { return HTAtomic128::Iterator{m_cells}; }

    Iterator end() const;

    size_t incrementGlobalOccupancy(size_t increment);

    size_t globalOccupancy() const;

    size_t capacity() const;

    Cells const &cells() const { return m_cells; }

    // Filled whenever occupancy reaches greater equal than half of the
    // capacity (alpha = 0.5 * capacity).
    bool filled() const;

    // Roams source to target with respect to the calling thread which will be given
    // a defined space of the source table to insert in target. This will guarantee
    // correctness and roams elements without any synchronisation between threads.
    // See chapter 5.3.2 for more details and proof.
    //
    // scale_factor: is power of 2
    //
    // p_count: must be >= 1
    void roam(ParallelHashMap::HTAtomic128 &target, uint32_t const p_count = 1,
              uint32_t const p_id = 0);

    std::pair<ParallelHashMap::HTAtomic128::Cells::const_iterator,
              ParallelHashMap::HTAtomic128::Cells::const_iterator>
    clusterRange(uint32_t const p_count = 1, uint32_t const p_id = 0);

private:
    static Cell makeCell(uint64_t key, uint64_t value) {
        Cell c;
        c.key = key;
        c.value = value;
        return c;
    }

    static Cell invalidCell() {
        return makeCell(ParallelHashMap::ht_invalid_key, ParallelHashMap::ht_invalid_value);
    }

    // Double Width Compare And Swap:
    //
    // If r matches expected, write desired to r and return true. Otherwise,
    // read r into expected and return false.
    //
    // We using this user defined dwcas since we can't be sure we're operating
    // atomically on two 64-bit values: https://timur.audio/dwcas-in-c
    static bool dwcas(HTAtomic128::Cell &r, HTAtomic128::Cell &expected,
                      HTAtomic128::Cell const &desired) {
        bool res;
#if __x86_64__
        asm volatile("lock cmpxchg16b %1"
                     : "=@ccz"(res), "+m"(r), "+a"(expected.key), "+d"(expected.value)
                     : "b"(desired.key), "c"(desired.value)
                     : "memory");
#else
        // For ARM64 we use the ld/stxp instruction pair to perform the DWCAS operation.
        // See:
        // https://github.com/boostorg/atomic/blob/master/include/boost/atomic/detail/core_arch_ops_gcc_aarch64.hpp
        HTAtomic128::Cell original;
        asm volatile("mov %w[success], #0\n\t"
                     "ldxp %x[original_0], %x[original_1], %[storage]\n\t"
                     "cmp %x[original_0], %x[expected_0]\n\t"
                     "ccmp %x[original_1], %x[expected_1], #0, eq\n\t"
                     "b.ne 1f\n\t"
                     "stxp %w[success], %x[desired_0], %x[desired_1], %[storage]\n\t"
                     "eor %w[success], %w[success], #1\n\t"
                     "1:\n\t"
                     : [success] "=&r"(res), [storage] "+Q"(r), [original_0] "=&r"(original.key),
                       [original_1] "=&r"(original.value)
                     : [desired_0] "r"(desired.key), [desired_1] "r"(desired.value),
                       [expected_0] "r"(expected.key), [expected_1] "r"(expected.value)
                     : "cc"
                       "memory");
#endif
        return res;
    }

    void moveCells(ParallelHashMap::HTAtomic128::Cells::const_iterator begin,
                   ParallelHashMap::HTAtomic128::Cells::const_iterator end,
                   ParallelHashMap::HTAtomic128 &target);

private:
    size_t m_capacity{ParallelHashMap::ht_begin_capacity};
    uint64_t m_log_capacity{static_cast<uint64_t>(std::log2(m_capacity))};
    Cells m_cells;
    std::atomic_size_t m_global_occupancy{0};
};

class ParallelHashMap::HTSyncData {
public:
    HTSyncData(std::unique_ptr<HTAtomic128> &_source, std::unique_ptr<HTAtomic128> &_target,
               std::atomic_uint32_t &_busy, std::atomic_bool &_request_growth, int _p_count,
               int _p_id, uint32_t _insert_threshold);

    std::unique_ptr<HTAtomic128> &source;
    std::unique_ptr<HTAtomic128> &target;

    // Represents a bitset where the bit is set to 1 at position p_id if the
    // thread is still busy operating on hash table operations.
    std::atomic_uint32_t &busy;
    uint32_t insert_counter{0};
    uint32_t const insert_threshold;

    // Whenever this is set, the hashtable has to be grown by all threads of the
    // omp team.
    std::atomic_bool &request_growth;
    int const p_count;
    int const p_id;
};

class ParallelHashMap::HTHandle {
public:
    HTHandle(HTSyncData sync_data);

    ~HTHandle();

    bool insert(uint64_t const key, uint64_t const value);

    uint64_t find(uint64_t const key) const;

    bool update(uint64_t const key, uint64_t const value) const;

    HTAtomic128 const &hashtable();

private:
    HTAtomic128 *m_ht;
    HTSyncData m_sync_data;
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_PARALLEL_HASH_MAP_HPP_
