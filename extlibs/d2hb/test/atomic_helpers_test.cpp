#include <catch2/catch_test_macros.hpp>

#include <d2hb/atomic_helpers.h>

TEST_CASE("atomic_helpers") {
    SECTION("unset_bit_atomically") {
        std::atomic_uint32_t fourth_bit_set{0b1111};

        d2hb::unset_bit_atomically(fourth_bit_set, 3);

        CHECK(fourth_bit_set == 0b0111);
    }

    SECTION("set_bit_atomically") {
        std::atomic_uint32_t fourth_bit_unset{0b0111};

        d2hb::set_bit_atomically(fourth_bit_unset, 3);

        CHECK(fourth_bit_unset == 0b1111);
    }

    SECTION("swap_atomics_non_atomically") {
        std::atomic_uint32_t a{0b0111};
        std::atomic_uint32_t a_copy{a.load()};
        std::atomic_uint32_t b{0b0110};
        std::atomic_uint32_t b_copy{b.load()};

        REQUIRE(a != b);
        REQUIRE(a == a_copy);
        REQUIRE(b == b_copy);

        d2hb::swap_atomics_non_atomically(a, b);

        CHECK(a == b_copy);
        CHECK(b == a_copy);
    }

    SECTION("increment atomically") {
        std::atomic_uint32_t a{1};
        uint32_t a_incremented_expected = a + 1;
        uint32_t a_incremented = d2hb::increment_atomically(a);

        CHECK(a_incremented_expected == a_incremented);
        CHECK(a_incremented_expected == a.load());
    }
}
