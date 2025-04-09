#pragma once

#include <d2hb/batcher.h>
#include <d2hb/block.h>
#include <d2hb/buckets.h>
#include <d2hb/dynamic_hashtable.h>
#include <d2hb/graph.h>

#include <omp.h>

#include <cassert>
#include <functional>
#include <memory>
#include <vector>

namespace d2hb {

template <typename E> class HypersparseMatrix {
  public:
    struct NeighborView {
        using iterator = typename BlockState<E>::iterator;
        using proxy = typename BlockState<E>::proxy;

        NeighborView(HypersparseMatrix* const m, size_t const idx) : m_g{m}, m_v_index{idx} {}

        auto begin() { return m_g->m_vertices[m_v_index].begin(); }

        auto end() { return m_g->m_vertices[m_v_index].valid_end(); }

        auto iterator_to(Vertex v) { return m_g->m_vertices[m_v_index].iterator_to(v); }

        bool exists(Vertex v) { return iterator_to(v) != end(); }

        void clear() { m_g->m_vertices[m_v_index].clear(); }

        size_t degree() { return end() - begin(); }

        proxy operator[](size_t i) { return *(begin() + i); }

        std::tuple<iterator, bool> insert(Vertex v, E ed) {
            auto& state = m_g->m_vertices[m_v_index];

            if (!state.full()) {
                return state.insert(v, ed);
            }

            // We need to reallocate the adjacency block.
            auto new_bhandle = m_g->m_manager->allocate_block(state.bsize() + 1);
            BlockState<E> new_block{new_bhandle, state};
            auto result = new_block.insert(v, ed);

            auto old_block = std::move(m_g->m_vertices[m_v_index]);
            auto old_bhandle = m_g->m_handles[m_v_index];
            m_g->m_vertices[m_v_index] = std::move(new_block);
            m_g->m_handles[m_v_index] = new_bhandle;
            m_g->m_manager->free_block(old_bhandle);

            return result;
        }

      private:
        HypersparseMatrix* const m_g;
        size_t const m_v_index;
    };

    struct ConstNeighborView {
        using iterator = typename BlockState<E>::const_iterator;
        using proxy = typename BlockState<E>::const_proxy;

        ConstNeighborView(HypersparseMatrix const* const g, size_t const idx)
            : m_g{g}, m_v_index{idx} {}

        auto begin() { return m_g->m_vertices[m_v_index].begin(); }

        auto end() { return m_g->m_vertices[m_v_index].valid_end(); }

        auto iterator_to(Vertex v) { return m_g->m_vertices[m_v_index].iterator_to(v); }

        bool exists(Vertex v) { return iterator_to(v) != end(); }

        size_t degree() { return end() - begin(); }

        proxy operator[](size_t i) { return *(begin() + i); }

      private:
        HypersparseMatrix const* const m_g;
        size_t const m_v_index;
    };

    friend void swap(HypersparseMatrix& p, HypersparseMatrix& q) noexcept {
        using std::swap;
        swap(p.m_vertices, q.m_vertices);
        swap(p.m_handles, q.m_handles);
        swap(p.m_manager, q.m_manager);
        swap(p.m_ht_custodian, q.m_ht_custodian);
    }

    HypersparseMatrix()
        : HypersparseMatrix(std::make_shared<BlockManager>(bytes_per_entry_for<E>())) {}

    HypersparseMatrix(std::shared_ptr<BlockManager> manager) : m_manager(std::move(manager)) {
        // Vertex at 0 is the "empty vertex"
        m_vertices.resize(1);
        m_handles.resize(1);
    }

    template <typename Edges_t> HypersparseMatrix(Edges_t&& edges) : HypersparseMatrix{} {
        auto handle = m_ht_custodian.make_handle();
        for (auto e = std::begin(edges); e != std::end(edges); ++e) {
            size_t stored_idx = handle->find(e->source);

            if (stored_idx == ht_invalid_value) {
                stored_idx = allocate_new_block(m_vertices, m_handles, 1);
                bool const inserted = handle->insert(e->source, stored_idx);
                assert(inserted);
            }

            auto result = neighbors_by_index(stored_idx).insert(e->target.vertex, e->target.data);
            assert(std::get<1>(result));
        }
    }

    HypersparseMatrix(const HypersparseMatrix& other) : HypersparseMatrix{} {
        if (other.vertices_count()) {
            resize(other.vertices_count());
        }

        m_ht_custodian = other.m_custodian;
    }

    HypersparseMatrix(HypersparseMatrix&& other) noexcept : HypersparseMatrix{} {
        swap(*this, other);
    }

    ~HypersparseMatrix() {
        for (auto bhandle : m_handles) {
            m_manager->free_block(bhandle);
        }
    }

    HypersparseMatrix& operator=(HypersparseMatrix other) {
        swap(*this, other);
        return *this;
    }

    void preallocate(Vertex u, Vertex degree) {
        auto handle = m_ht_custodian.make_handle();

        size_t stored_idx = handle->find(u);

        if (stored_idx == ht_invalid_value) {
            stored_idx = allocate_new_block(m_vertices, m_handles, degree);
            handle->insert(u, stored_idx);
        } else {
            auto& associated_block = m_vertices[stored_idx];

            if (associated_block.bsize() < degree) {
                auto new_bhandle = m_manager->allocate_block(degree);
                BlockState<E> new_block{new_bhandle, associated_block};

                auto old_block = std::move(m_vertices[u]);
                auto old_bhandle = m_handles[stored_idx];
                m_vertices[stored_idx] = std::move(new_block);
                m_handles[stored_idx] = new_bhandle;
                m_manager->free_block(old_bhandle);
            }
        }
    }

    bool insert(Vertex u, Vertex v, E ed, bool do_update = true) {
        auto handle = m_ht_custodian.make_handle();
        size_t stored_idx = handle->find(u);

        if (stored_idx == ht_invalid_value) {
            stored_idx = allocate_new_block(m_vertices, m_handles, 1);
            bool const inserted = handle->insert(u, stored_idx);
            assert(inserted);
        }

        auto result = neighbors_by_index(stored_idx).insert(v, ed);
        if (std::get<1>(result)) {
            return true;
        }

        if (do_update) {
            std::get<0>(result)->data() = ed;
            return true;
        }

        return false;
    }

    // Use this method to apply concurrent operations while inserting some local
    // (begin, end) batch of edges. This will reduce the synchronisation
    // complexity.
    //
    // Using this method, there shall be no edge source index u handled by any
    // other thread but this. Otherwise, the call of this method yields
    // undefined behaviour.
    template <typename EdgeIt, typename DataFn, typename OnUpdate>
    void insert(EdgeIt begin, EdgeIt end, DataFn&& data_fn, bool update, OnUpdate on_update) {
        constexpr uint64_t highest_bit = uint64_t(1) << 63;
        auto is_local_idx = [](uint64_t const idx) -> bool { return idx & highest_bit; };

        auto ht_handle = m_ht_custodian.make_handle();

        std::vector<BlockState<E>> local_vertices;
        std::vector<BlockHandle> local_handles;
        std::vector<Vertex> local_ids;

        for (auto edge = begin; edge != end; ++edge) {
            uint64_t stored_idx = ht_handle->find(edge->source);

            if (stored_idx == ht_invalid_value) {
                stored_idx = allocate_new_block(local_vertices, local_handles, 1);
                stored_idx = stored_idx | highest_bit;
                bool const inserted = ht_handle->insert(edge->source, stored_idx);
                assert(inserted);
                local_ids.push_back(edge->source);
            }

            bool const is_local = is_local_idx(stored_idx);
            if (is_local) {
                stored_idx = stored_idx & ~highest_bit;
            }

            std::vector<BlockState<E>>& vertices = is_local ? local_vertices : m_vertices;
            std::vector<BlockHandle>& handles = is_local ? local_handles : m_handles;

            BlockState<E>& state = vertices[stored_idx];
            if (!state.full()) {
                auto result = state.insert(edge->target.vertex, data_fn(edge));
                if (update && !std::get<1>(result)) {
                    on_update(std::get<0>(result)->data(), data_fn(edge));
                }
            } else {
                // We need to reallocate the adjacency block.
                auto new_bhandle = m_manager->allocate_block(state.bsize() + 1);
                BlockState<E> new_block{new_bhandle, state};
                auto result = new_block.insert(edge->target.vertex, data_fn(edge));
                if (update && !std::get<1>(result)) {
                    on_update(std::get<0>(result)->data(), data_fn(edge));
                }

                auto old_block = std::move(vertices[stored_idx]);
                auto old_bhandle = handles[stored_idx];
                vertices[stored_idx] = std::move(new_block);
                handles[stored_idx] = new_bhandle;
                m_manager->free_block(old_bhandle);
            }
        }

#pragma omp barrier

        uint64_t current_end = 0;

#pragma omp critical
        {
            current_end = m_vertices.size();
            // Move vertex states and handles to the end of the global
            // containers
            m_vertices.insert(m_vertices.end(), std::make_move_iterator(local_vertices.begin()),
                              std::make_move_iterator(local_vertices.end()));
            m_handles.insert(m_handles.end(), std::make_move_iterator(local_handles.begin()),
                             std::make_move_iterator(local_handles.end()));
        }

        // Now update all vertex IDs with their new index within the
        // containers
        for (auto id : local_ids) {
            bool const updated = ht_handle->update(id, current_end);
            assert(updated);
            ++current_end;
        }
    }

    // Use this function for concurrent map operations on a batch of edges. The
    // apply_fn will be provided the existant neighbour (vertex with |N(v)| > 0
    // - non zero row), and the current edge iterator.
    template <typename EdgeIt, typename ApplyFn>
    void map(EdgeIt begin, EdgeIt end, ApplyFn apply_fn) {
        auto table = m_ht_custodian.current_table();

        for (auto edge = begin; edge != end; ++edge) {
            size_t stored_idx = table->find(edge->source);

            if (stored_idx != ht_invalid_value) {
                apply_fn(neighbors_by_index(stored_idx), edge);
            }
        }
    }

    bool removeEdge(Vertex u, Vertex target) {
        size_t const idx = m_ht_custodian.current_table()->find(u);

        if (idx == ht_invalid_value) {
            return false;
        }

        auto& associated_block = m_vertices[idx];
        if (associated_block.remove(target)) {
            return true;
        }

        return false;
    }

    template <typename L> void sort(Vertex u, L less) {
        size_t const idx = m_ht_custodian.current_table()->find(u);

        if (idx != ht_invalid_value) {
            m_vertices[idx].sort(std::move(less));
        }
    }

    NeighborView neighbors(Vertex const u) {
        auto handle = m_ht_custodian.make_handle();
        auto idx = handle->find(u);

        if (idx != ht_invalid_value) {
            return neighbors_by_index(idx);
        }

        size_t const stored_idx = allocate_new_block(m_vertices, m_handles, 1);
        handle->insert(u, stored_idx);
        return neighbors_by_index(stored_idx);
    }

    ConstNeighborView neighbors(Vertex const u) const {
        auto ht = m_ht_custodian.current_table();

        size_t const idx = ht->find(u);

        if (idx != ht_invalid_value) {
            return neighbors_by_index(idx);
        }

        // This index is reserved for all vertices with |n(u)| = 0
        return neighbors_by_index(0);
    }

    template <typename F> void for_nodes(F&& handle) const {
        for (auto entry : *m_ht_custodian.current_table()) {
            handle(entry.key);
        }
    }

    template <typename F> void for_edges(F&& handle) const {
        for (auto entry : *m_ht_custodian.current_table()) {
            if (m_vertices[entry.value].degree() > 0) {
                auto& block = m_vertices[entry.value];
                for (auto it = block.cbegin(); it != block.valid_end(); it++) {
                    handle(entry.key, it->vertex(), it->data());
                }
            }
        }
    }

    template <typename F> void for_neighbors_node(Vertex const u, F&& handle) const {
        size_t const idx = m_ht_custodian.current_table()->find(u);

        if (idx == ht_invalid_value) {
            return;
        }

        auto& block = m_vertices[idx];

        for (auto it = block.cbegin(); it != block.valid_end(); it++) {
            handle(it->vertex());
        }
    }

    template <typename F> void for_neighbors_target(Vertex const u, F&& handle) const {
        size_t const idx = m_ht_custodian.current_table()->find(u);

        if (idx == ht_invalid_value) {
            return;
        }

        auto& block = m_vertices[idx];

        for (auto it = block.cbegin(); it != block.valid_end(); it++) {
            handle(it->vertex(), it->data());
        }
    }

    Vertex degree(Vertex u) const {
        size_t const idx = m_ht_custodian.current_table()->find(u);

        if (idx == ht_invalid_value) {
            return Vertex(0);
        }

        return m_vertices[idx].degree();
    }

    // Returns the maximum known vertex ID + 1. Adding + 1 since we assume that
    // the vertex IDs start from 0.
    Vertex vertices_count() const {
        auto const table = m_ht_custodian.current_table();

        Vertex max_v = 0;
        for (auto c : *table) {
            if (c.key > max_v) {
                max_v = c.key;
            }
        }

        return max_v + 1;
    }

    Vertex edges_count() const {
        return std::accumulate(
            std::cbegin(m_vertices), std::cend(m_vertices), 0u,
            [](Vertex sum, const BlockState<E>& block) { return sum + block.degree(); });
    }

  private:
    NeighborView neighbors_by_index(size_t idx) { return {this, idx}; }

    ConstNeighborView neighbors_by_index(size_t idx) const { return {this, idx}; }

    size_t allocate_new_block(std::vector<BlockState<E>>& vertices,
                              std::vector<BlockHandle>& handles, Vertex degree) {
        d2hb::BlockHandle bhandle = m_manager->allocate_block(degree);
        BlockState<E> block{bhandle};

        vertices.emplace_back(std::move(block));
        handles.emplace_back(bhandle);

        return vertices.size() - 1;
    }

  private:
    HTCustodian m_ht_custodian;
    std::vector<BlockState<E>> m_vertices;
    std::vector<BlockHandle> m_handles;
    std::shared_ptr<BlockManager> m_manager;
}; // namespace d2hb

template <typename E, class Hashmap, class HashmapCell>
inline std::tuple<Degree, Vertex> max_degree(const HypersparseMatrix<E>& matrix) {
    Vertex v = 0;
    for (Vertex u = 1; u < matrix.vertices_count(); ++u) {
        if (matrix.degree(u) > matrix.degree(v)) {
            v = u;
        }
    }
    assert(v < matrix.vertices_count());
    return {matrix.degree(v), v};
}

} // namespace d2hb
