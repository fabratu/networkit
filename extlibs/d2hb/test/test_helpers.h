#pragma once

#include <d2hb/dynamic_hashtable.h>
#include <d2hb/graph.h>

#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <algorithm>

using MockupData = std::vector<std::pair<d2hb::Vertex, size_t>>;

inline MockupData generate_mockup_data(size_t const fill_64bit_values = 148) {
    MockupData mockup_data;

    constexpr d2hb::Vertex minimum_of_range = 1;
    constexpr d2hb::Vertex maximum_of_range = 356;
    auto random_range =
        Catch::Generators::RandomIntegerGenerator<d2hb::Vertex>(minimum_of_range, maximum_of_range);

    d2hb::Vertex u = 0;
    for (size_t i = 0; i < fill_64bit_values; ++i) {
        random_range.next();
        u += random_range.get();
        mockup_data.push_back(std::make_pair(u, i));
    }

    return mockup_data;
}

inline bool check_entries_for_mockup_data(MockupData const& mockup_data,
                                          d2hb::HTAtomic128 const& ht) {
    bool all_good = true;
    for (size_t i = 0; i < mockup_data.size() && all_good; ++i) {
        size_t const retrieved_idx = ht.find(mockup_data[i].first);

        all_good = retrieved_idx < mockup_data.size();
        all_good = all_good && retrieved_idx == mockup_data[i].second;
    }
    return all_good;
}

inline bool get_and_check_all(std::vector<d2hb::Vertex> const& mockup_data, d2hb::HTAtomic128& ht) {
    size_t idx = 0;
    bool all_good = true;
    for (size_t i = 0; i < mockup_data.size() && all_good; ++i) {
        idx = ht.find(mockup_data[i]);

        all_good = idx < mockup_data.size();
        all_good = all_good && mockup_data[idx] == mockup_data[i];
    }

    return all_good;
}

inline bool insert_data(d2hb::HTAtomic128& ht, MockupData const& data) {
    bool all_inserted_successful = true;
    for (size_t i = 0; i < data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = ht.insert(data[i].first, data[i].second);
    }
    return all_inserted_successful;
}

inline bool insert_data(d2hb::HTHandle& handle, MockupData const& data) {
    bool all_inserted_successful = true;
    for (size_t i = 0; i < data.size() && all_inserted_successful; ++i) {
        all_inserted_successful = handle.insert(data[i].first, data[i].second);
    }
    return all_inserted_successful;
}
