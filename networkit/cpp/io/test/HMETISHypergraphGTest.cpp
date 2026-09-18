/*
 * HMETISHypergraphGTest.cpp
 */

#include <gtest/gtest.h>

#include <string>
#include <unordered_set>

#include <networkit/graph/Hypergraph.hpp>
#include <networkit/io/HMETISHypergraphReader.hpp>
#include <networkit/io/HMETISHypergraphWriter.hpp>

namespace NetworKit {

class HMETISHypergraphGTest : public testing::Test {};

namespace {

void expectEqual(const Hypergraph &actual, const Hypergraph &expected) {
    ASSERT_EQ(actual.numberOfNodes(), expected.numberOfNodes());
    ASSERT_EQ(actual.numberOfEdges(), expected.numberOfEdges());
    ASSERT_EQ(actual.isWeighted(), expected.isWeighted());
    ASSERT_EQ(actual.isIncidenceWeighted(), expected.isIncidenceWeighted());

    expected.forEdges([&](edgeid eid) {
        ASSERT_TRUE(actual.hasEdge(eid));
        EXPECT_EQ(actual.edgeMembers(eid), expected.edgeMembers(eid));
        EXPECT_DOUBLE_EQ(actual.getEdgeWeight(eid), expected.getEdgeWeight(eid));
        for (node u : expected.edgeMembers(eid)) {
            EXPECT_DOUBLE_EQ(actual.getIncidenceWeight(u, eid),
                             expected.getIncidenceWeight(u, eid));
        }
    });
}

} // namespace

TEST_F(HMETISHypergraphGTest, testReadAugmentedExample) {
    const Hypergraph hypergraph =
        HMETISHypergraphReader{}.read("input/edge-bar-reviews-swapped.hmetis");

    EXPECT_EQ(hypergraph.numberOfNodes(), 15);
    EXPECT_EQ(hypergraph.numberOfEdges(), 222);
    EXPECT_TRUE(hypergraph.isWeighted());
    EXPECT_TRUE(hypergraph.isIncidenceWeighted());

    EXPECT_DOUBLE_EQ(hypergraph.getEdgeWeight(0), 1.0);
    EXPECT_EQ(hypergraph.edgeMembers(0), (std::unordered_set<node>{0, 1}));
    EXPECT_DOUBLE_EQ(hypergraph.getIncidenceWeight(0, 0), 23.0);
    EXPECT_DOUBLE_EQ(hypergraph.getIncidenceWeight(1, 0), 11.0);
    EXPECT_DOUBLE_EQ(hypergraph.weightedDegree(0), 5572.0);

    edgeweight edgeWeightSum = 0.0;
    edgeweight incidenceWeightSum = 0.0;
    count numberOfIncidences = 0;
    hypergraph.forEdges([&](edgeid eid, edgeweight weight) {
        edgeWeightSum += weight;
        for (node u : hypergraph.edgeMembers(eid)) {
            incidenceWeightSum += hypergraph.getIncidenceWeight(u, eid);
            ++numberOfIncidences;
        }
    });
    EXPECT_DOUBLE_EQ(edgeWeightSum, 1234.0);
    EXPECT_DOUBLE_EQ(incidenceWeightSum, 3313.0);
    EXPECT_EQ(numberOfIncidences, 346);
}

TEST_F(HMETISHypergraphGTest, testRoundTripAllSupportedFormats) {
    for (bool edgeWeighted : {false, true}) {
        for (bool incidenceWeighted : {false, true}) {
            Hypergraph expected{5, 0, edgeWeighted, incidenceWeighted};
            const edgeid first = expected.addEdge({0, 2, 4});
            const edgeid second = expected.addEdge({1, 2});
            expected.addEdge({});

            if (edgeWeighted) {
                expected.setEdgeWeight(first, 2.5);
                expected.setEdgeWeight(second, -3.25);
            }
            if (incidenceWeighted) {
                expected.setIncidenceWeight(0, first, 4.5);
                expected.setIncidenceWeight(2, first, 0.25);
                expected.setIncidenceWeight(4, first, -2.0);
                expected.setIncidenceWeight(1, second, 8.0);
                expected.setIncidenceWeight(2, second, 1.5);
            }

            const std::string path =
                "output/hmetis-roundtrip-" + std::to_string(static_cast<int>(edgeWeighted)) + "-"
                + std::to_string(static_cast<int>(incidenceWeighted)) + ".hmetis";
            HMETISHypergraphWriter{}.write(expected, path);
            expectEqual(HMETISHypergraphReader{}.read(path), expected);
        }
    }
}

TEST_F(HMETISHypergraphGTest, testExampleRoundTrip) {
    const Hypergraph expected =
        HMETISHypergraphReader{}.read("input/edge-bar-reviews-swapped.hmetis");
    const std::string path = "output/edge-bar-reviews-swapped.hmetis";

    HMETISHypergraphWriter{}.write(expected, path);
    expectEqual(HMETISHypergraphReader{}.read(path), expected);
}

} // namespace NetworKit
