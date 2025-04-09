#include <catch2/catch_test_macros.hpp>

#include <d2hb/submatrix.h>

using namespace d2hb;
using ED = d2hb::EdgeData;
using Tar = d2hb::Target;

class TestMatrix {
  public:
    TestMatrix() {
        // Matrix A 3x3:
        // [  1 | 3 | 0 ]
        // [  4 | 5 | 8 ]
        // [ 13 | 4 | 0 ]
        // clang-format off
        d2hb::Edges edges_a = { 
            {0, Tar{0, ED{1.0}}},  {0, Tar{1, ED{3.0}}}, {0, Tar{2, ED{0.0}}},
            {1, Tar{0, ED{4.0}}},  {1, Tar{1, ED{5.0}}}, {1, Tar{2, ED{8.0}}},
            {2, Tar{0, ED{13.0}}}, {2, Tar{1, ED{4.0}}}, {2, Tar{2, ED{0.0}}},
        };
        // clang-format on

        a = d2hb::Matrix<d2hb::EdgeData>(std::move(edges_a));
    }

  protected:
    d2hb::Matrix<d2hb::EdgeData> a;
};

TEST_CASE_METHOD(TestMatrix, "construct with matrix", "[construction]") {
    d2hb::Submatrix<d2hb::EdgeData> a_sub(a);

    CHECK(a_sub.edges_count() == 0);
    CHECK(a_sub.set().size() == 0);
    CHECK(&a_sub.matrix() == &a);
}

TEST_CASE_METHOD(TestMatrix, "construct with matrix, begin, and end vertex", "[construction]") {
    d2hb::Submatrix<d2hb::EdgeData> a_sub(a, 0, a.vertices_count());

    CHECK(a_sub.edges_count() == 9);
    CHECK(a_sub.set().size() == 3);
    CHECK(&a_sub.matrix() == &a);
}

TEST_CASE_METHOD(TestMatrix, "extend irregular", "[extend_set]") {
    d2hb::Submatrix<d2hb::EdgeData> a_sub(a, 0, 1);
    CHECK(a_sub.set()[0] == 0);
    CHECK(a_sub.edges_count() == 3);

    a_sub.extend_set(2, 3);

    CHECK(a_sub.edges_count() == 6);
    CHECK(a_sub.set()[0] == 0);
    CHECK(a_sub.set()[1] == 2);
}

TEST_CASE_METHOD(TestMatrix, "construct with SubmatrixSet", "[construction]") {
    d2hb::SubmatrixSet set{0, 1};
    Vertex const set_size = set.size();

    d2hb::Submatrix<d2hb::EdgeData> a_sub(a, std::move(set));

    CHECK(a_sub.edges_count() == 6);
    CHECK(a_sub.set().size() == set_size);
    CHECK(&a_sub.matrix() == &a);
}

TEST_CASE_METHOD(TestMatrix, "copy constructor", "[construction]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a(a, 0, a.vertices_count());
    d2hb::Submatrix<d2hb::EdgeData> sub_a_copy(sub_a);

    CHECK(sub_a.set() == sub_a_copy.set());
    CHECK(&sub_a.matrix() == &sub_a_copy.matrix());
}

TEST_CASE_METHOD(TestMatrix, "move constructor", "[construction]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a_moved(a, 0, a.vertices_count());
    d2hb::Submatrix<d2hb::EdgeData> sub_a_copy(sub_a_moved);
    d2hb::Submatrix<d2hb::EdgeData> sub_a(std::move(sub_a_moved));

    CHECK(&sub_a.matrix() == &sub_a_copy.matrix());
    CHECK(sub_a.set() == sub_a_copy.set());
}

TEST_CASE_METHOD(TestMatrix, "copy assignment constructor", "[construction]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a(a, 0, a.vertices_count());
    d2hb::Submatrix<d2hb::EdgeData> sub_a_copy = sub_a;

    CHECK(sub_a.set() == sub_a_copy.set());
    CHECK(&sub_a.matrix() == &sub_a_copy.matrix());
}

TEST_CASE_METHOD(TestMatrix, "move assignment constructor", "[construction]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a_moved(a, 0, a.vertices_count());
    d2hb::Submatrix<d2hb::EdgeData> sub_a_copy(sub_a_moved);
    d2hb::Submatrix<d2hb::EdgeData> sub_a = std::move(sub_a_moved);

    CHECK(sub_a.set() == sub_a_copy.set());
    CHECK(&sub_a.matrix() == &sub_a_copy.matrix());
}

TEST_CASE_METHOD(TestMatrix, "operator==", "[comparison]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a(a, 0, a.vertices_count());
    d2hb::Submatrix<d2hb::EdgeData> sub_a_copy = sub_a;

    CHECK(sub_a == sub_a_copy);
}

TEST_CASE_METHOD(TestMatrix, "operator!=", "[comparison]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a(a, 0, a.vertices_count());
    REQUIRE(sub_a.vertices_count() == 3);
    d2hb::Submatrix<d2hb::EdgeData> sub_b(a, 0, 2);
    REQUIRE(sub_b.vertices_count() == 2);

    CHECK(sub_a != sub_b);
}

TEST_CASE_METHOD(TestMatrix, "for_neighbors_node()", "[for_neighbors_node]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a(a, 0, a.vertices_count());
    std::vector<bool> active_nodes(3, false);

    sub_a.for_neighbors_node(0, [&](d2hb::Vertex const v) { active_nodes[v] = true; });

    CHECK(active_nodes[0]);
    CHECK(active_nodes[1]);
    CHECK(active_nodes[2]);
}

TEST_CASE_METHOD(TestMatrix, "neighbors(u)", "[neighbors]") {
    d2hb::Submatrix<d2hb::EdgeData> sub_a(a, 0, a.vertices_count());

    CHECK(sub_a.neighbors(0).iterator_to(2)->data().weight ==
          a.neighbors(0).iterator_to(2)->data().weight);
}

TEST_CASE_METHOD(TestMatrix, "vector of submatrices", "[std::vector]") {
    std::vector<d2hb::Submatrix<d2hb::EdgeData>> submatrices(2, a);
    REQUIRE(submatrices.size() == 2);

    CHECK(submatrices[0].edges_count() == 0);
    CHECK(submatrices[1].edges_count() == 0);

    CHECK(submatrices[0].set().size() == 0);
    CHECK(submatrices[1].set().size() == 0);
}

TEST_CASE_METHOD(TestMatrix, "vector of submatrices using function initialisation",
                 "[lambda][init]") {
    std::vector<d2hb::Submatrix<d2hb::EdgeData>> submatrices = [&]() {
        std::vector<d2hb::Submatrix<d2hb::EdgeData>> submatrices(2, a);
        submatrices[0].extend_set(0, 2);
        submatrices[1].extend_set(1, 3);
        return submatrices;
    }();

    REQUIRE(submatrices.size() == 2);

    CHECK(submatrices[0].edges_count() == 6);
    CHECK(submatrices[1].edges_count() == 6);

    CHECK(submatrices[0].set().size() == 2);
    CHECK(submatrices[1].set().size() == 2);
}
