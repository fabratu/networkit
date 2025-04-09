#include <iostream>

#include <d2hb/batcher.h>
#include <d2hb/dynamic_doubly_hashed_blocks.h>
#include <d2hb/graph.h>

#include <omp.h>

#include <limits>
#include <random>

using DataPack = std::vector<std::pair<d2hb::Vertex, size_t>>;
using namespace d2hb;
using HyperGraph = HypersparseMatrix<EdgeData>;

DataPack generate_data(std::mt19937& generator, size_t const fill_count = 148) {
    DataPack data;

    constexpr d2hb::Vertex minimum_of_range = 1;
    constexpr d2hb::Vertex maximum_of_range = 356;
    std::uniform_int_distribution<uint32_t> m_dist{minimum_of_range, maximum_of_range};

    d2hb::Vertex u = 0;
    for (size_t i = 0; i < fill_count; ++i) {
        u += m_dist(generator);
        data.push_back(std::make_pair(u, i));
    }

    return data;
}

std::chrono::milliseconds ht_routine(std::mt19937& generator, unsigned int const thread_count) {
    constexpr size_t fill_count = 1'000'000;
    DataPack data = generate_data(generator, fill_count);

    d2hb::HTCustodian custodian{32u};

    auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel num_threads(thread_count) shared(custodian)
    {
        uint32_t const my_id = omp_get_thread_num();
        auto handle = custodian.make_handle();
        auto it = std::cbegin(data) + my_id;
        while (it != std::cend(data)) {
            if (handle->find(it->first) == ht_invalid_value) {
                handle->insert(it->first, it->second);
            }

            std::advance(it, std::min(thread_count, uint32_t(std::distance(it, std::cend(data)))));
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Filled cells: " << custodian.current_table()->global_occupancy() << std::endl;

    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
}

Edges generate_edges(size_t const fill_count, std::mt19937& generator) {
    Edges edges;
    constexpr Vertex minimum_of_range = 0;
    constexpr Vertex maximum_of_range = ht_invalid_value - 1;
    std::uniform_int_distribution<Vertex> m_dist{minimum_of_range, maximum_of_range};
    std::uniform_real_distribution<Weight> m_float_dist{1.0, 10.0};

    EdgeID id = 0;
    for (size_t i = 0; i < fill_count; ++i) {
        edges.push_back(Edge{m_dist(generator),
                             Target{m_dist(generator), EdgeData{m_float_dist(generator), id}}});
        id++;
    }

    return edges;
}

bool all_edges_exist(HyperGraph const& graph, Edges const& edges) {
    bool all_edges_exist = true;
    for (auto e = std::cbegin(edges); e != std::cend(edges) && all_edges_exist; ++e) {
        all_edges_exist = all_edges_exist && graph.neighbors(e->source).exists(e->target.vertex);
    }
    return all_edges_exist;
}

std::chrono::milliseconds gds_routine(std::mt19937& generator, unsigned int const thread_count) {
    constexpr size_t fill_count = 1'000'000;

    Edges edges = generate_edges(fill_count, generator);

    auto cmp = [](Edge const& a, Edge const& b) { return a.source < b.source; };
    std::sort(std::begin(edges), std::end(edges), std::move(cmp));

    HyperGraph graph;

    HypersparseBatchParallelizer<Edges, Edges::const_iterator> par;
    auto start = std::chrono::high_resolution_clock::now();

    auto batch_begin = std::cbegin(edges);
    auto batch_end = std::cbegin(edges);
    size_t const batch_size = 100'000;
    while (batch_begin != std::cend(edges)) {
        // We advance batch end but only until we reached the end
        std::advance(batch_end,
                     std::min(batch_size, size_t(std::distance(batch_end, std::cend(edges)))));

        auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
        auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
            auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
            graph.insert(begin, end, std::move(data_fn), false, std::move(on_update));
        };

        par.apply_exclusive(batch_begin, batch_end, std::move(fun));
        batch_begin = batch_end;
    }

    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Edge count: " << graph.edges_count() << std::endl;

    std::cout << "All edges exist: " << all_edges_exist(graph, edges) << std::endl;

    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
}

std::chrono::milliseconds small_example(std::mt19937& generator, unsigned int const thread_count) {
    Edges first_batch;
    first_batch.push_back(Edge{0, Target{1, EdgeData{1.0f, 0}}});
    first_batch.push_back(Edge{0, Target{2, EdgeData{1.0f, 1}}});
    first_batch.push_back(Edge{1, Target{3, EdgeData{1.0f, 2}}});
    first_batch.push_back(Edge{1, Target{4, EdgeData{1.0f, 3}}});

    Edges second_batch;
    second_batch.push_back(Edge{0, Target{4, EdgeData{1.0f, 4}}});
    second_batch.push_back(Edge{0, Target{5, EdgeData{1.0f, 5}}});
    second_batch.push_back(Edge{2, Target{3, EdgeData{1.0f, 6}}});

    second_batch.push_back(Edge{3, Target{1, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{2, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{3, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{4, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{5, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{6, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{7, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{8, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{9, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{10, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{11, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{12, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{13, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{14, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{15, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{16, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{17, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{18, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{19, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{20, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{21, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{22, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{23, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{24, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{25, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{26, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{27, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{28, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{29, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{30, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{31, EdgeData{1.0f, 7}}});
    second_batch.push_back(Edge{3, Target{32, EdgeData{1.0f, 7}}});

    using HyperGraph = HypersparseMatrix<EdgeData>;
    HyperGraph graph;

    auto on_update = [](EdgeData& current, EdgeData&& update) { current = update; };
    auto fun = [&](Edges::const_iterator begin, Edges::const_iterator end) {
        auto data_fn = [](Edges::const_iterator edge) { return edge->target.data; };
        graph.insert(begin, end, std::move(data_fn), false, std::move(on_update));
    };

    HypersparseBatchParallelizer<Edges, Edges::const_iterator> par;

    auto start = std::chrono::high_resolution_clock::now();
    par.apply_exclusive(std::begin(first_batch), std::end(first_batch), std::move(fun));

    std::cout << "Inserting second BATCH" << std::endl;

    par.apply_exclusive(std::begin(second_batch), std::end(second_batch), std::move(fun));
    auto stop = std::chrono::high_resolution_clock::now();

    std::cout << "Edge count: " << graph.edges_count() << std::endl;
    std::cout << "All edges exist: "
              << (all_edges_exist(graph, first_batch) && all_edges_exist(graph, second_batch))
              << std::endl;

    return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
}

int main(int argc, char** argv) {
    unsigned int const thread_count = 4;
    omp_set_num_threads(thread_count);

    std::mt19937 generator;

    auto duration = gds_routine(generator, thread_count);

    std::cout << "Duration: " << duration.count() << std::endl;

    return 0;
}
