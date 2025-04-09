#include <catch2/catch_test_macros.hpp>

#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include "graph_io.h"
#include "test_helpers.h"

#include <d2hb/dynamic_hashed_blocks.h>
#include <d2hb/dynamic_hashtable.h>

#include <vector>

using namespace d2hb;

// https://networkrepository.com/bio-celegans.php
static std::string graph_name{"tango.edges"};
static std::string graph_path =
#ifdef D2HB_TEST_GRAPH_DIR
    std::string(D2HB_TEST_GRAPH_DIR) + "/"
#else
    "test/graphs/"
#endif
    + graph_name;

TEST_CASE("fast_mod") {
    uint32_t idx = 8;
    size_t capacity = 16;

    uint32_t idx_wrapped = fast_mod(idx, capacity);

    CHECK(idx_wrapped == idx);

    idx = 13;
    idx_wrapped = fast_mod(idx, capacity);

    CHECK(idx_wrapped == idx);

    idx = 17;
    idx_wrapped = fast_mod(idx, capacity);

    CHECK(idx_wrapped == 1);
}

TEST_CASE("key_space") { CHECK(ht_key_space == 64); }

TEST_CASE("Using hashtable with Hypersparse Tango Graph") {
    Edges edges = read_graph_unweighted(graph_path, false);

    // + 2 = 18 but the last two vertices can not be counted since vertex 16,
    //   17, and 18 are not part of the edges.
    unsigned int const tango_vertex_count = 16;
    unsigned int const tango_edge_count = 14;

    Degree const tango_max_degree = 2;
    std::vector<Degree> extracted_degrees = degrees_from(edges);

    REQUIRE(edges.size() == tango_edge_count);
    REQUIRE(graph::vertex_count(edges) == tango_vertex_count);

    // Representing a vector data structure holding all targets in consecutive
    // order just like the dhb matrix class.
    std::vector<Target> targets;

    std::for_each(std::begin(edges), std::end(edges), [&](Edge const& e) {
        targets.push_back(
            Target{e.target.vertex, EdgeData{e.target.data.weight, e.target.data.id}});
    });

    REQUIRE(targets.size() == edges.size());

    SECTION("Sanity Check") {
        HTAtomic128 ht{ht_begin_capacity};

        ht.insert(targets[0].vertex, 0);
        CHECK(ht.find(targets[0].vertex) == 0);
    }

    SECTION("Manual Test") {
        HTAtomic128 ht{ht_begin_capacity};

        size_t idx = 0;
        d2hb::Vertex source = 5;
        d2hb::Vertex target = 6;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        size_t target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 6;
        target = 5;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 7;
        target = 8;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 8;
        target = 7;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 9;
        target = 10;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 10;
        target = 9;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 11;
        target = 12;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 12;
        target = 11;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 13;
        target = 14;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 14;
        target = 13;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 13;
        target = 15;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        CHECK(!ht.insert(source, idx));

        idx++;
        source = 15;
        target = 13;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        ht.insert(source, idx);
        target_index = ht.find(source);
        CHECK(targets[target_index].vertex == target);

        idx++;
        source = 14;
        target = 15;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        CHECK(!ht.insert(source, idx));

        idx++;
        source = 15;
        target = 14;
        REQUIRE(edges[idx].source == source);
        REQUIRE(targets[idx].vertex == target);
        CHECK(!ht.insert(source, idx));

        idx++;
        REQUIRE(idx == edges.size());

        // What if we want to get a not yet added source
        CHECK(ht.find(18) == ht_invalid_value);
    }

    SECTION("Automatically Insert and Get All") {
        HTAtomic128 ht{ht_begin_capacity};

        d2hb::Vertex source = 0;
        d2hb::Vertex target = 0;
        size_t target_index = 0;

        bool already_inserted_are_valid = true;
        bool find_succeeded = true;

        for (size_t i = 0; i < edges.size() && already_inserted_are_valid && find_succeeded; ++i) {
            source = edges[i].source;
            target = targets[i].vertex;

            if (!ht.insert(source, i)) {
                already_inserted_are_valid = ht.find(source) != ht_invalid_value;
            } else {
                target_index = ht.find(source);

                find_succeeded = edges[target_index].source == source;
                find_succeeded = find_succeeded && targets[target_index].vertex == target;
            }
        }
    }

    SECTION("update") {
        HTAtomic128 ht;

        d2hb::Vertex source = 0;
        d2hb::Vertex target = 0;
        size_t target_index = 0;

        bool already_inserted_are_valid = true;
        bool find_succeeded = true;

        for (size_t i = 0; i < edges.size() && already_inserted_are_valid && find_succeeded; ++i) {
            source = edges[i].source;
            target = targets[i].vertex;

            if (!ht.insert(source, i)) {
                already_inserted_are_valid = ht.find(source) != ht_invalid_value;
            } else {
                target_index = ht.find(source);

                find_succeeded = edges[target_index].source == source;
                find_succeeded = find_succeeded && targets[target_index].vertex == target;
            }
        }

        bool all_updated_successfully = true;
        for (size_t i = 0; i < edges.size() && all_updated_successfully; ++i) {
            all_updated_successfully =
                all_updated_successfully && ht.update(edges[i].source, targets[i].vertex + 1);
        }
        CHECK(all_updated_successfully);

        for (size_t i = 0; i < edges.size() && all_updated_successfully; ++i) {
            all_updated_successfully = ht.find(source) == targets[i].vertex + 1;
        }
    }
}

TEST_CASE("generate_mockup_data") {
    constexpr size_t size = 6844;
    MockupData const mockup_data = generate_mockup_data(size);

    REQUIRE(mockup_data.size() == size);
}

TEST_CASE("Hashfolklore, Capacity") {
    SECTION("1 MB") {
        MockupData const mockup_data = generate_mockup_data(125000);

        HTAtomic128 ht{131072};

        for (size_t i = 0; i < mockup_data.size(); ++i) {
            ht.insert(mockup_data[i].first, mockup_data[i].second);
        }

        REQUIRE(check_entries_for_mockup_data(mockup_data, ht));
    }
}

TEST_CASE("HTAtomic128, Misc") {
    SECTION("Empty") {
        HTAtomic128 ht{ht_begin_capacity};

        CHECK(ht.find(0) == ht_invalid_value);
        CHECK(ht.find(155) == ht_invalid_value);
    }
}

TEST_CASE("HTAtomic128, [Copy, Move, Assign, etc.]") {
    std::vector<d2hb::Vertex> mockup_data;
    constexpr size_t fill_value_with_64bit_values = 148;
    constexpr d2hb::Vertex vertex_offset_from_index = 25;

    for (size_t i = 0; i < fill_value_with_64bit_values; ++i) {
        mockup_data.push_back(i + vertex_offset_from_index);
    }

    HTAtomic128 ht{ht_begin_capacity};

    for (size_t i = 0; i < mockup_data.size(); ++i) {
        ht.insert(mockup_data[i], i);
    }

    SECTION("Copy Constructor") {
        HTAtomic128 dvh_copy(ht);

        bool all_good = get_and_check_all(mockup_data, dvh_copy);

        REQUIRE(get_and_check_all(mockup_data, ht));
        CHECK(get_and_check_all(mockup_data, ht));
    }

    SECTION("Move Constructor") {
        HTAtomic128 dvh_moved(std::move(ht));

        REQUIRE(get_and_check_all(mockup_data, dvh_moved));
    }

    SECTION("swap") {
        std::vector<d2hb::Vertex> mockup_data_b;
        constexpr size_t fill_value_with_64bit_values_b = 54;
        constexpr d2hb::Vertex vertex_offset_from_index_b = 26;

        for (size_t i = 0; i < fill_value_with_64bit_values_b; ++i) {
            mockup_data_b.push_back(i + vertex_offset_from_index_b);
        }

        HTAtomic128 ht_b{ht_begin_capacity};

        for (size_t i = 0; i < mockup_data_b.size(); ++i) {
            ht_b.insert(mockup_data_b[i], i);
        }

        ht_b.increment_global_occupancy(mockup_data_b.size());

        swap(ht, ht_b);

        CHECK(get_and_check_all(mockup_data_b, ht));
        CHECK(get_and_check_all(mockup_data, ht_b));
    }

    SECTION("Copy assignment") {
        HTAtomic128 ht_b{ht_begin_capacity};
        ht_b = ht;

        CHECK(get_and_check_all(mockup_data, ht_b));
        CHECK(get_and_check_all(mockup_data, ht));
    }
}

TEST_CASE("HTAtomic128, Input Iterator") {
    MockupData const mockup_data = generate_mockup_data();

    HTAtomic128 const ht = [&]() {
        HTAtomic128 ht{ht_begin_capacity};
        bool all_inserted_successful = true;
        for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
            all_inserted_successful = ht.insert(mockup_data[i].first, mockup_data[i].second);
        }

        REQUIRE(all_inserted_successful);

        return ht;
    }();

    REQUIRE(check_entries_for_mockup_data(mockup_data, ht));

    SECTION("begin, end") {
        HTAtomic128::Iterator it = ht.begin();
        HTAtomic128::Iterator it_end = ht.end();

        CHECK(it != it_end);

        auto it_2 = ht.begin();

        CHECK(it == it_2);

        auto it_end_2 = ht.end();

        CHECK(it_end == it_end_2);

        it.invalidate();

        CHECK(it == it_end);
    }

    SECTION("manual check") {
        HTAtomic128::Iterator it = ht.begin();
        size_t iterated_elements = 0;

        bool all_valid = true;
        bool not_end = true;
        while (iterated_elements < mockup_data.size() && all_valid && not_end) {
            all_valid = it->key != ht_invalid_key;
            not_end = it != ht.end();
            it++;
            iterated_elements++;
        }

        CHECK(all_valid);
        CHECK(not_end);
        CHECK(iterated_elements == mockup_data.size());

        it++;
        CHECK(it == ht.end());
    }

    SECTION("iterate") {
        bool keys_are_correct = true;
        bool values_are_correct = true;
        size_t elements_iterated = 0;

        bool it_not_end = true;
        for (auto it = ht.begin(); it != ht.end() && keys_are_correct && values_are_correct; ++it) {
            keys_are_correct = it->key == mockup_data[it->value].first;
            values_are_correct = it->value == mockup_data[it->value].second;
            it_not_end = it != ht.end();
            elements_iterated++;
        }

        CHECK(it_not_end);
        CHECK(keys_are_correct);
        CHECK(values_are_correct);
        CHECK(elements_iterated == mockup_data.size());
    }

    SECTION("for") {
        size_t elements_iterated = 0;

        bool keys_are_correct = true;
        bool values_are_correct = true;
        for (auto c : ht) {
            keys_are_correct = keys_are_correct && (c.key == mockup_data[c.value].first);
            values_are_correct = values_are_correct && (c.value == mockup_data[c.value].second);
            elements_iterated++;
        }

        CHECK(keys_are_correct);
        CHECK(values_are_correct);
        CHECK(elements_iterated == mockup_data.size());
    }
}

TEST_CASE("HTAtomic128, Concurrent access") {
    unsigned int const thread_count = 2;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    SECTION("1 MB") {
        constexpr size_t one_megabyte_filled_with_64bit_values = 125000;

        auto mockup_data = generate_mockup_data(one_megabyte_filled_with_64bit_values);

        constexpr size_t begin_capacity = 131072;
        HTAtomic128 ht{begin_capacity};

#pragma omp parallel for
        for (auto e : mockup_data) {
            ht.insert(e.first, e.second);
        }

        CHECK(check_entries_for_mockup_data(mockup_data, ht));
    }

    SECTION("Random writes") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HTAtomic128 ht;

#pragma omp parallel
        {
            int thread_id = omp_get_thread_num();

            if (thread_id == 0) {
                ht.insert(0, 0);
                ht.insert(1, 0);
            }

            if (thread_id == 1) {
                ht.insert(0, 1);
                ht.insert(1, 1);
            }
        }

        Vertex zero = ht.find(0);
        CHECK((zero == 0 || zero == 1));
        Vertex one = ht.find(1);
        CHECK((one == 0 || one == 1));
    }

    SECTION("Linearizeabilty") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HTAtomic128 ht;

#pragma omp parallel
        {
            int thread_id = omp_get_thread_num();

            if (thread_id == 0) {
                ht.insert(0, 0);
                ht.insert(1, 0);
            }

#pragma omp barrier
            if (thread_id == 1) {
                CHECK(ht.find(0) == 0);
                CHECK(ht.find(1) == 0);
            }
        }
    }

    SECTION("Update") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HTAtomic128 ht;

        ht.insert(0, 0);
        ht.insert(1, 0);

#pragma omp parallel
        {
            int thread_id = omp_get_thread_num();

            if (thread_id == 0) {
                ht.update(0, 2);
            } else {
                ht.update(0, 3);
            }
        }

        uint64_t const updated_value = ht.find(0);

        CHECK((updated_value == 2 || updated_value == 3));
    }
}

TEST_CASE("HTAtomic128, insert, dual thread") {
    unsigned int const thread_count = 2;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generate_mockup_data(fill_64bit_values);

    constexpr size_t begin_capacity = 16384u;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = true;
#pragma omp parallel
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = ht_first.insert(mockup_data[i].first, mockup_data[i].second);
    }

    CHECK(all_inserted_successful);
    CHECK(check_entries_for_mockup_data(mockup_data, ht_first));
}

TEST_CASE("hashtable capacity greater alpha") {
    size_t capacity = 32;
    size_t occupancy = 10;

    CHECK(!ht_filled(occupancy, capacity));

    occupancy = 15;
    CHECK(!ht_filled(occupancy, capacity));

    occupancy = 16;
    CHECK(!ht_filled(occupancy, capacity));

    occupancy = 17;
    CHECK(ht_filled(occupancy, capacity));

    occupancy = capacity * 2;
    CHECK(ht_filled(occupancy, capacity));
}
