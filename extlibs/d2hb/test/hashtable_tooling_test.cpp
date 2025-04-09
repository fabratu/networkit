#include <catch2/catch_test_macros.hpp>

#include <d2hb/dynamic_hashtable.h>
#include <d2hb/graph.h>

#include "test_helpers.h"

#include <iostream>
#include <vector>

using namespace d2hb;

void print(HTAtomic128& ht) {
    for (auto c : ht.cells()) {

        auto pair = [](HTAtomic128::Cell const& c) {
            if (c.key == ht_invalid_key) {
                return std::make_pair(std::string("I"), std::string("V"));
            }

            return std::make_pair(std::to_string(c.key), std::to_string(c.value));
        }(c);

        std::cout << "[" << pair.first << "," << pair.second << "]"
                  << "\n";
    }

    std::cout << std::endl;
}

TEST_CASE("HTAtomic128", "cluster_range") {
    constexpr size_t fill_64bit_values = 345;
    MockupData mockup_data = generate_mockup_data(fill_64bit_values);

    constexpr size_t begin_capacity = 512;
    HTAtomic128 source{begin_capacity};

    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = source.insert(mockup_data[i].first, mockup_data[i].second);
    }
    REQUIRE(all_inserted_successful);
    REQUIRE(check_entries_for_mockup_data(mockup_data, source));

    auto const cells_begin = std::cbegin(source.cells());
    auto const cells_end = std::cend(source.cells());

    SECTION("for single thread") {
        auto c_range = cluster_range(cells_begin, cells_end);
        CHECK(c_range.first == cells_begin);
        CHECK(c_range.second == cells_end);
    }

    SECTION("for dual thread") {
        auto c_range_t0 = cluster_range(cells_begin, cells_end, 2, 0);
        CHECK(c_range_t0.first == cells_begin);
        CHECK(c_range_t0.second != cells_end);

        auto c_range_t1 = cluster_range(cells_begin, cells_end, 2, 1);
        CHECK(c_range_t1.first != cells_begin);
        CHECK(c_range_t1.second == cells_end);

        CHECK(c_range_t0.second == c_range_t1.first);
    }

    SECTION("for quadruple thread") {
        uint32_t const thread_count = 4;

        auto c_range_t0 = cluster_range(cells_begin, cells_end, thread_count, 0);
        CHECK(c_range_t0.first == cells_begin);
        CHECK(c_range_t0.second != cells_end);

        auto c_range_t1 = cluster_range(cells_begin, cells_end, thread_count, 1);
        CHECK(c_range_t1.first != cells_begin);
        CHECK(c_range_t1.second != cells_end);
        CHECK(c_range_t0.second == c_range_t1.first);

        auto c_range_t2 = cluster_range(cells_begin, cells_end, thread_count, 2);
        CHECK(c_range_t2.first != cells_begin);
        CHECK(c_range_t2.second != cells_end);
        CHECK(c_range_t1.second == c_range_t2.first);

        auto c_range_t3 = cluster_range(cells_begin, cells_end, thread_count, 3);
        CHECK(c_range_t3.first != cells_begin);
        CHECK(c_range_t3.second == cells_end);
        CHECK(c_range_t2.second == c_range_t3.first);
    }
}

TEST_CASE("HTAtomic128, Roaming, single thread") {
    constexpr size_t fill_64bit_values = 23;
    MockupData mockup_data = generate_mockup_data(fill_64bit_values);

    constexpr size_t begin_capacity = 32;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = true;
    for (size_t i = 0; i < mockup_data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = ht_first.insert(mockup_data[i].first, mockup_data[i].second);
    }

    REQUIRE(all_inserted_successful);

    REQUIRE(check_entries_for_mockup_data(mockup_data, ht_first));

    constexpr size_t next_capacity = begin_capacity << 1;
    HTAtomic128 ht_second{next_capacity};

    roam(ht_first, ht_second);

    CHECK(check_entries_for_mockup_data(mockup_data, ht_second));
}

TEST_CASE("HTAtomic128, Roaming, dual thread") {
    unsigned int const thread_count = 2;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    constexpr size_t fill_64bit_values = 23;
    MockupData mockup_data = generate_mockup_data(fill_64bit_values);

    constexpr size_t begin_capacity = 32;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = insert_data(ht_first, mockup_data);

    REQUIRE(all_inserted_successful);
    REQUIRE(check_entries_for_mockup_data(mockup_data, ht_first));

    constexpr size_t next_capacity = begin_capacity << 1;
    HTAtomic128 ht_second{next_capacity};

#pragma omp parallel
    { roam(ht_first, ht_second, thread_count, omp_get_thread_num()); }

    CHECK(check_entries_for_mockup_data(mockup_data, ht_second));
}

TEST_CASE("HTAtomic128, Big Roaming, four thread") {
    unsigned int const thread_count = 4;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    constexpr size_t fill_64bit_values = 10254;
    MockupData mockup_data = generate_mockup_data(fill_64bit_values);

    constexpr size_t begin_capacity = 16384u;
    HTAtomic128 ht_first{begin_capacity};

    bool all_inserted_successful = insert_data(ht_first, mockup_data);

    REQUIRE(all_inserted_successful);
    REQUIRE(check_entries_for_mockup_data(mockup_data, ht_first));

    constexpr size_t next_capacity = begin_capacity << 1;
    HTAtomic128 ht_second{next_capacity};

#pragma omp parallel
    { roam(ht_first, ht_second, thread_count, omp_get_thread_num()); }

    CHECK(check_entries_for_mockup_data(mockup_data, ht_second));
}

TEST_CASE("Tooling, set_bit_atomically") {
    std::atomic_uint32_t bitmask{0u};

    uint32_t desired = 1u;

    uint32_t result = set_bit_atomically(bitmask, 0);

    CHECK(desired == result);
}

TEST_CASE("HTAtomic128, HTHandle, HTRoaming") {
    MockupData const mockup_data = generate_mockup_data();
    std::mt19937 generator;

    SECTION("insert all, single thread") {
        unsigned int const thread_count = 1;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        std::unique_ptr<HTAtomic128> source = std::make_unique<HTAtomic128>(ht_begin_capacity);
        std::unique_ptr<HTAtomic128> target;

        std::atomic_uint32_t busy_bitset{0u};
        std::atomic_bool request_growth{false};

        HTSyncData sync_data{
            source, target, busy_bitset, request_growth, thread_count, omp_get_thread_num(), 2};
        HTHandle handle{sync_data};

        bool const all_inserted_successful = insert_data(handle, mockup_data);

        REQUIRE(all_inserted_successful);

        CHECK(check_entries_for_mockup_data(mockup_data, handle.hashtable()));

        CHECK(handle.hashtable().global_occupancy() == mockup_data.size());
    }

    SECTION("Using two threads to insert.") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        std::unique_ptr<HTAtomic128> source = std::make_unique<HTAtomic128>(ht_begin_capacity);
        std::unique_ptr<HTAtomic128> target;

        std::atomic_uint32_t busy_bitset{0u};
        std::atomic_bool request_growth{false};

#pragma omp parallel
        {
            HTSyncData sync_data{
                source, target, busy_bitset, request_growth, thread_count, omp_get_thread_num(), 2};
            HTHandle handle{sync_data};

            auto mockup_data_begin = std::begin(mockup_data);
            std::advance(mockup_data_begin, omp_get_thread_num() == 0 ? 0 : mockup_data.size() / 2);

            auto mockup_data_end = std::begin(mockup_data);
            std::advance(mockup_data_end,
                         omp_get_thread_num() == 0 ? mockup_data.size() / 2 : mockup_data.size());

            bool successful_inserts = true;
            for (auto it = mockup_data_begin; it != mockup_data_end && successful_inserts; ++it) {
                bool success = handle.insert(it->first, it->second);
            }

            CHECK(successful_inserts);
        }

        HTAtomic128& ht = *source.get();

        CHECK(check_entries_for_mockup_data(mockup_data, ht));

        CHECK(ht.global_occupancy() == mockup_data.size());
    }

    SECTION("Dual Thread Insert, Roaming Required") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        std::unique_ptr<HTAtomic128> source = std::make_unique<HTAtomic128>(32u);
        std::unique_ptr<HTAtomic128> target;

        size_t const mockup_data_share = mockup_data.size() / thread_count;

        std::atomic_uint32_t busy_bitset{0u};
        std::atomic_bool request_growth{false};

#pragma omp parallel
        {
            int const thread_id = omp_get_thread_num();
            HTSyncData sync_data{source,    target, busy_bitset, request_growth, thread_count,
                                 thread_id, 2};
            HTHandle handle{sync_data};

            auto begin = std::begin(mockup_data);
            std::advance(begin, mockup_data_share * thread_id);

            auto end = begin;

            if (thread_id == (thread_count - 1)) {
                end = std::end(mockup_data);
            } else {
                std::advance(end, mockup_data_share);
            }

            bool all_inserted_successful = true;
            for (auto it = begin; it != end && all_inserted_successful; it++) {
                all_inserted_successful = handle.insert(it->first, it->second);
            }

            CHECK(all_inserted_successful);
        }

        HTAtomic128& ht = *source.get();

        CHECK(check_entries_for_mockup_data(mockup_data, ht));
    }
}

TEST_CASE("HTCustodian, Copy, Assign, Move, etc") {
    unsigned int const thread_count = 1;
    omp_set_num_threads(thread_count);
    REQUIRE(omp_get_max_threads() == thread_count);

    SECTION("Swap") {

        HTCustodian custodian_a{32};
        HTCustodian custodian_b{32};

        MockupData mockup_data = generate_mockup_data();

        auto handle_a = custodian_a.make_handle();
        auto handle_b = custodian_b.make_handle();

        REQUIRE(handle_a->insert(mockup_data[0].first, mockup_data[0].second));
        REQUIRE(handle_b->insert(mockup_data[1].first, mockup_data[1].second));

        swap(custodian_a, custodian_b);

        auto handle_a_new = custodian_a.make_handle();
        auto handle_b_new = custodian_b.make_handle();

        CHECK(handle_a_new->find(mockup_data[1].first) == mockup_data[1].second);
        CHECK(handle_b_new->find(mockup_data[0].first) == mockup_data[0].second);
    }

    SECTION("Copy Constructor") {
        HTCustodian custodian_a;
        auto mockup_data = generate_mockup_data();

        auto handle_a = custodian_a.make_handle();

        bool all_inserted_successful = insert_data(*handle_a.get(), mockup_data);

        CHECK(check_entries_for_mockup_data(mockup_data, handle_a->hashtable()));

        HTCustodian custodian_b(custodian_a);

        auto handle_b = custodian_b.make_handle();

        CHECK(check_entries_for_mockup_data(mockup_data, handle_b->hashtable()));
    }

    SECTION("Move Constructor") {
        HTCustodian custodian_a;
        auto mockup_data = generate_mockup_data();

        auto handle_a = custodian_a.make_handle();

        bool all_inserted_successful = insert_data(*handle_a.get(), mockup_data);

        HTCustodian custodian_b(std::move(custodian_a));

        auto handle_b = custodian_b.make_handle();

        CHECK(check_entries_for_mockup_data(mockup_data, handle_b->hashtable()));
    }

    SECTION("Copy assignment") {
        HTCustodian custodian_a;
        auto mockup_data = generate_mockup_data();

        auto handle_a = custodian_a.make_handle();

        bool all_inserted_successful = insert_data(*handle_a.get(), mockup_data);

        HTCustodian custodian_b;
        custodian_b = custodian_a;

        auto handle_b = custodian_b.make_handle();

        CHECK(check_entries_for_mockup_data(mockup_data, handle_b->hashtable()));
    }
}

TEST_CASE("HTAtomic128, HTCustodian, Threading") {
    SECTION("Single Thread") {
        unsigned int const thread_count = 1;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HTCustodian custodian{32};

        constexpr size_t fill_64bit_values = 10254;
        MockupData mockup_data = generate_mockup_data(fill_64bit_values);

#pragma omp parallel num_threads(thread_count) shared(custodian)
        {
            auto handle = custodian.make_handle();
            for (auto it = std::cbegin(mockup_data); it != std::cend(mockup_data); ++it) {
                handle->insert(it->first, it->second);
            }
        }

        CHECK(check_entries_for_mockup_data(mockup_data, *custodian.current_table()));
    }

    SECTION("Dual Thread") {
        unsigned int const thread_count = 2;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HTCustodian custodian{32};

        constexpr size_t fill_64bit_values = 10254;
        MockupData mockup_data = generate_mockup_data(fill_64bit_values);

#pragma omp parallel num_threads(thread_count) shared(custodian, mockup_data)
        {
            uint32_t const thread_id = omp_get_thread_num();
            auto handle = custodian.make_handle();

            auto local_begin = std::cbegin(mockup_data);
            size_t const elements_per_thread = mockup_data.size() / thread_count;

            if (thread_id != 0) {
                std::advance(local_begin, thread_id * elements_per_thread);
            }

            auto local_end = local_begin;

            if (thread_id == (thread_count - 1)) {
                local_end = std::cend(mockup_data);
            } else {
                std::advance(local_end, (thread_id + 1) * elements_per_thread);
            }

            bool all_inserted = true;
            for (auto it = local_begin; it != local_end && all_inserted; ++it) {
                all_inserted = handle->insert(it->first, it->second);
            }

            CHECK(all_inserted);
        }

        CHECK(check_entries_for_mockup_data(mockup_data, *custodian.current_table()));
    }

    SECTION("Quadruple Thread") {
        unsigned int const thread_count = 4;
        omp_set_num_threads(thread_count);
        REQUIRE(omp_get_max_threads() == thread_count);

        HTCustodian custodian{32};

        constexpr size_t fill_64bit_values = 10254;
        MockupData mockup_data = generate_mockup_data(fill_64bit_values);

#pragma omp parallel num_threads(thread_count) shared(custodian, mockup_data)
        {
            uint32_t const thread_id = omp_get_thread_num();
            auto handle = custodian.make_handle();

            auto local_begin = std::cbegin(mockup_data);
            size_t const elements_per_thread = mockup_data.size() / thread_count;

            if (thread_id != 0) {
                std::advance(local_begin, thread_id * elements_per_thread);
            }

            auto local_end = local_begin;

            if (thread_id == (thread_count - 1)) {
                local_end = std::cend(mockup_data);
            } else {
                std::advance(local_end, elements_per_thread);
            }

            bool all_inserted = true;
            size_t counter = 0;
            for (auto it = local_begin; it != local_end; ++it) {
                all_inserted = all_inserted && handle->insert(it->first, it->second);
                if (!all_inserted) {
                    counter++;
                }
            }

            CHECK(counter == 0);

            CHECK(all_inserted);
        }

        CHECK(check_entries_for_mockup_data(mockup_data, *custodian.current_table()));
    }
}
