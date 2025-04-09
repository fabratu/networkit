#include <catch2/catch_test_macros.hpp>

#include "graph_io.h"

#include <d2hb/dynamic_doubly_hashed_blocks.h>
#include <d2hb/graph.h>

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <experimental/filesystem>
#include <map>

namespace fs = std::experimental::filesystem;

using namespace d2hb;

// https://networkrepository.com/bio-celegans.php
static std::string graph_name{"bio-celegans.mtx"};
static std::string graph_path =
#ifdef D2HB_TEST_GRAPH_DIR
    std::string(D2HB_TEST_GRAPH_DIR) + "/"
#else
    "test/graphs/"
#endif
    + graph_name;

// we start at 0 so we must have 453 + 1 vertices and not 453 as described on
// https://networkrepository.com/bio-celegans.php
constexpr Vertex bio_celegans_vertex_count = 453 + 1;

TEST_CASE("HypersparseMatrix, basic functionality") {
    unsigned int const thread_count = 1;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    Edges edges = read_graph_unweighted(graph_path, false);
    unsigned int const bio_celegans_edge_count = edges.size();

    Degree const bio_celegans_max_degree = 237;
    std::vector<Degree> extracted_degrees = degrees_from(edges);

    SECTION("sanity check") {
        HypersparseMatrix<EdgeData> m(std::move(edges));
        CHECK(m.vertices_count() == bio_celegans_vertex_count);

        size_t const max_degree_vertex = std::distance(
            std::begin(extracted_degrees),
            std::max_element(std::begin(extracted_degrees), std::end(extracted_degrees)));

        CHECK(m.edges_count() == bio_celegans_edge_count);
        CHECK(m.degree(max_degree_vertex) == bio_celegans_max_degree);
    }

    SECTION("NeighborView") {
        HypersparseMatrix<EdgeData> m(std::move(edges));
        // N(89) = 13, 31
        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 35);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        auto n89 = nv.begin();
        CHECK(n89 != nv.end());

        // yes, we assume an order here
        CHECK(n89->vertex() == 13);
        CHECK(++n89 != nv.end());
        CHECK(n89->vertex() == 31);
        CHECK(++n89 != nv.end());
        CHECK(n89->vertex() == 86);

        auto n31 = nv.iterator_to(31);
        CHECK(n31->vertex() == 31);

        CHECK(nv[0].vertex() == 13);
        std::tuple<HypersparseMatrix<EdgeData>::NeighborView::iterator, bool> r =
            nv.insert(14, EdgeData{10.f, 450});

        CHECK(std::get<0>(r)->vertex() == 14);
        CHECK(nv.degree() == 36);
        CHECK(nv.exists(14));
    }

    SECTION("ConstNeighborView") {
        HypersparseMatrix<EdgeData> const m(std::move(edges));
        // N(89) = 13, 31
        HypersparseMatrix<EdgeData>::ConstNeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == 35);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        auto n89 = nv.begin();
        CHECK(n89 != nv.end());

        // yes, we assume an order here
        CHECK(n89->vertex() == 13);
        REQUIRE(++n89 != nv.end());
        CHECK(n89->vertex() == 31);
        REQUIRE(++n89 != nv.end());
        CHECK(n89->vertex() == 86);

        auto n31 = nv.iterator_to(31);
        CHECK(n31->vertex() == 31);

        CHECK(nv[0].vertex() == 13);
    }

    SECTION("ConstNeighborView", "non existent vertex") {
        HypersparseMatrix<EdgeData> const m(std::move(edges));

        HypersparseMatrix<EdgeData>::ConstNeighborView nv = m.neighbors(8042);
        CHECK(nv.degree() == 0);
        CHECK(!nv.exists(13));
        CHECK(!nv.exists(31));
        CHECK(!nv.exists(86));

        auto n89 = nv.begin();
        CHECK(n89 == nv.end());

        auto n31 = nv.iterator_to(31);
        CHECK(n31 == nv.end());
    }

    SECTION("for_nodes") {
        std::vector<Vertex> vertices_from_edges = [&]() {
            std::vector<Vertex> vertices;
            for (auto e : edges) {
                vertices.push_back(e.source);
            }

            std::sort(std::begin(vertices), std::end(vertices));
            auto last = std::unique(std::begin(vertices), std::end(vertices));
            vertices.erase(last, std::end(vertices));

            return vertices;
        }();

        HypersparseMatrix<EdgeData> const m(std::move(edges));

        std::vector<Vertex> vertices;

        m.for_nodes([&](Vertex u) { vertices.push_back(u); });

        std::sort(std::begin(vertices), std::end(vertices));

        REQUIRE(vertices.size() == vertices_from_edges.size());

        bool nodes_are_equivalent = true;
        for (size_t i = 0; i < vertices.size() && nodes_are_equivalent; ++i) {
            nodes_are_equivalent = vertices[i] == vertices_from_edges[i];
        }

        CHECK(nodes_are_equivalent);
    }

    SECTION("for_edges") {
        std::map<std::tuple<Vertex, Vertex>, bool> edge_map = [&]() {
            std::map<std::tuple<Vertex, Vertex>, bool> edge_map;
            for (auto e : edges) {
                auto u_v = std::make_tuple(e.source, e.target.vertex);
                edge_map.insert({u_v, false});
            }

            return edge_map;
        }();

        HypersparseMatrix<EdgeData> const m(std::move(edges));

        m.for_edges([&](Vertex u, Vertex v, EdgeData data) {
            auto u_v = std::make_tuple(u, v);
            edge_map[u_v] = true;
        });

        bool all_edges_visited = true;

        for (auto e : edge_map) {
            all_edges_visited = all_edges_visited && e.second;
        }

        CHECK(all_edges_visited);
    }

    SECTION("sort") {
        HypersparseMatrix<EdgeData> m(std::move(edges));
        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);

        Vertex v89_degree = 35;

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        nv.insert(e89_14.target.vertex, e89_14.target.data);
        v89_degree++;
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};
        nv.insert(e89_8.target.vertex, e89_8.target.data);
        v89_degree++;

        m.sort(89, [](BlockState<EdgeData>::Entry& a, BlockState<EdgeData>::Entry& b) {
            return a.vertex < b.vertex;
        });

        CHECK(nv.degree() == v89_degree);

        auto n = nv.begin();
        CHECK(n->vertex() == 8);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 13);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 14);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 31);
        REQUIRE(++n != nv.end());
        CHECK(n->vertex() == 86);

        // Iterate over the rest of all neighbours
        for (unsigned int i = 0; i < v89_degree - 5; ++i) {
            n++;
        }

        CHECK(++n == nv.end());
    }

    SECTION("removeEdge") {
        HypersparseMatrix<EdgeData> m(std::move(edges));
        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);

        Vertex v89_degree = 35;

        CHECK(nv.degree() == v89_degree);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        CHECK(m.removeEdge(89, 13));
        v89_degree--;

        CHECK(nv.degree() == v89_degree);

        CHECK(!nv.exists(13));

        CHECK(!m.removeEdge(89, 13));
    }
}

TEST_CASE("HypersparseMatrix, insert, update, single thread") {
    unsigned int const thread_count = 1;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    Edges edges = read_graph_unweighted(graph_path, false);

    SECTION("insert all, no update") {
        HypersparseMatrix<EdgeData> m;

        Vertex v89_degree = 35;

        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), false, std::move(on_update));
        };
        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.apply_non_exclusive(std::begin(edges), std::end(edges), std::move(fun));

        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == v89_degree);
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("insert partly, no update") {
        HypersparseMatrix<EdgeData> m(std::move(edges));

        Vertex v89_degree = 35;

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edges new_edges{e89_14, e89_8};

        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), false, std::move(on_update));
        };
        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.apply_non_exclusive(new_edges.begin(), new_edges.end(), std::move(fun));

        v89_degree += 2;

        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == v89_degree);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("insert with update, single thread") {
        HypersparseMatrix<EdgeData> m(std::move(edges));

        Vertex v89_degree = 35;
        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edge e89_13_update{89, Target{13, EdgeData{11.f, 455}}};

        Edges new_edges{e89_14, e89_8, e89_13_update};

        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto key = [](Edge e) { return e.source; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), true, std::move(on_update));
        };

        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.sort_and_apply_exclusive(new_edges.begin(), new_edges.end(), std::move(key),
                                     std::move(cmp), std::move(fun));

        v89_degree += 2;

        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == v89_degree);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        CHECK(nv.iterator_to(e89_13_update.target.vertex)->data().weight ==
              e89_13_update.target.data.weight);
    }
}

TEST_CASE("batcher") {
    Edges edges = read_graph_unweighted(graph_path, false);

    SECTION("non_exclusive_thread_batch, quad thread") {
        auto batch_0 = non_exclusive_thread_batch(std::cbegin(edges), std::cend(edges), 4, 0);
        auto batch_1 = non_exclusive_thread_batch(std::cbegin(edges), std::cend(edges), 4, 1);
        auto batch_2 = non_exclusive_thread_batch(std::cbegin(edges), std::cend(edges), 4, 2);
        auto batch_3 = non_exclusive_thread_batch(std::cbegin(edges), std::cend(edges), 4, 3);

        CHECK(std::get<0>(batch_0) == std::cbegin(edges));
        CHECK(std::get<1>(batch_0) == std::get<0>(batch_1));
        CHECK(std::get<1>(batch_1) == std::get<0>(batch_2));
        CHECK(std::get<1>(batch_2) == std::get<0>(batch_3));
        CHECK(std::get<1>(batch_3) == std::cend(edges));
    }
}

TEST_CASE("HypersparseMatrix, Concurrency") {
    SECTION("parallel insert no update dual thread") {
        Edges edges = read_graph_unweighted(graph_path, false);
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HypersparseMatrix<EdgeData> m(std::move(edges));

        Vertex v89_degree = 35;

        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edges new_edges{e89_14, e89_8};

        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto key = [](Edge e) { return e.source; };
        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), false, std::move(on_update));
        };

        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;

        par.sort_and_apply_exclusive(new_edges.begin(), new_edges.end(), std::move(key),
                                     std::move(cmp), std::move(fun));

        v89_degree += 2;

        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == v89_degree);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));
    }

    SECTION("parallel insert with update dual thread") {
        Edges edges = read_graph_unweighted(graph_path, false);

        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HypersparseMatrix<EdgeData> m(std::move(edges));

        Vertex v89_degree = 35;
        Edge e89_14{89, Target{14, EdgeData{10.f, 450}}};
        Edge e89_8{89, Target{8, EdgeData{10.f, 451}}};

        Edge e89_13_update{89, Target{13, EdgeData{11.f, 455}}};

        Edges new_edges{e89_14, e89_8, e89_13_update};

        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto key = [](Edge e) { return e.source; };
        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), true, std::move(on_update));
        };

        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.sort_and_apply_exclusive(new_edges.begin(), new_edges.end(), std::move(key),
                                     std::move(cmp), std::move(fun));

        v89_degree += 2;

        HypersparseMatrix<EdgeData>::NeighborView nv = m.neighbors(89);
        CHECK(nv.degree() == v89_degree);
        CHECK(nv.exists(14));
        CHECK(nv.exists(8));
        CHECK(nv.exists(13));
        CHECK(nv.exists(31));
        CHECK(nv.exists(86));

        CHECK(nv.iterator_to(e89_13_update.target.vertex)->data().weight ==
              e89_13_update.target.data.weight);
    }

    SECTION("parallel insert, using batch paralleliser, dual thread") {
        Edges edges = read_graph_unweighted(graph_path, false);

        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HypersparseMatrix<EdgeData> m;
        auto on_update = [](EdgeData& current, EdgeData&& update) { current = current; };
        auto key = [](Edge e) { return e.source; };
        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), false, std::move(on_update));
        };

        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.sort_and_apply_exclusive(std::begin(edges), std::end(edges), std::move(key),
                                     std::move(cmp), std::move(fun));

        CHECK(m.edges_count() == edges.size());
        CHECK(m.vertices_count() == bio_celegans_vertex_count);
    }

    SECTION("parallel insert, using batch paralleliser, quadruple thread") {
        Edges edges = read_graph_unweighted(graph_path, false);

        unsigned int const thread_count = 4;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HypersparseMatrix<EdgeData> m;
        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto key = [](Edge& e) { return e.source; };
        auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.insert(begin, end, std::move(data_fn), false, std::move(on_update));
        };

        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.sort_and_apply_exclusive(std::begin(edges), std::end(edges), std::move(key),
                                     std::move(cmp), std::move(fun));

        CHECK(m.edges_count() == edges.size());
        CHECK(m.vertices_count() == bio_celegans_vertex_count);
    }

    SECTION("parallel update, dual thread") {
        Edges edges = read_graph_unweighted(graph_path, false);

        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        Edges const update_edges = [&]() {
            Edges updates;
            for (auto e : edges) {
                updates.push_back({e.source, Target{e.target.vertex, {2.f, e.target.data.id}}});
            }
            return updates;
        }();

        HypersparseMatrix<EdgeData> m(std::move(edges));
        REQUIRE(m.edges_count() == edges.size());
        REQUIRE(m.vertices_count() == bio_celegans_vertex_count);

        auto apply_fn = [](HypersparseMatrix<EdgeData>::NeighborView&& nv,
                           Edges::const_iterator edge) {
            auto entry = nv.iterator_to(edge->target.vertex);

            if (entry != nv.end()) {
                entry->data().weight = edge->target.data.weight;
            }
        };

        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            m.map(begin, end, std::move(apply_fn));
        };

        HypersparseBatchParallelizer<Edge, Edges::const_iterator> par;
        par.apply_non_exclusive(std::begin(update_edges), std::end(update_edges), std::move(fun));

        bool all_data_correct = true;
        m.for_edges([&](Vertex u, Vertex v, EdgeData const& data) {
            all_data_correct = all_data_correct && data.weight == 2.f;
        });

        CHECK(all_data_correct);
    }
}
