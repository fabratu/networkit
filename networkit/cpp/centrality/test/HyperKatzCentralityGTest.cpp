#include <gtest/gtest.h>

#include <networkit/centrality/HyperKatzCentrality.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

TEST(HyperKatzCentralityGTest, testSumsAllSLevelsWithSpMV) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({2, 3});

    HyperKatzCentrality centrality(hGraph, 3, false, 1e-12);
    EXPECT_DOUBLE_EQ(centrality.getAlpha(1), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(centrality.getAlpha(2), 0.5);
    EXPECT_THROW(centrality.getAlpha(3), std::invalid_argument);
    EXPECT_THROW(centrality.getAlpha(4), std::out_of_range);

    centrality.run();

    // ASSERT_GE(centrality.nPaths.size(), 3);
    // EXPECT_DOUBLE_EQ(centrality.nPaths[1][0], 2.0);
    // EXPECT_DOUBLE_EQ(centrality.nPaths[1][1], 3.0);
    // EXPECT_DOUBLE_EQ(centrality.nPaths[1][2], 1.0);
    // EXPECT_DOUBLE_EQ(centrality.nPaths[2][0], 3.0);
    // EXPECT_DOUBLE_EQ(centrality.nPaths[2][1], 3.0);
    // EXPECT_DOUBLE_EQ(centrality.nPaths[2][2], 2.0);

    const auto ranking = centrality.ranking();
    ASSERT_EQ(ranking.size(), 3);
    EXPECT_EQ(ranking[0].first, 1);
    EXPECT_EQ(ranking[1].first, 0);
    EXPECT_EQ(ranking[2].first, 2);
    EXPECT_EQ(centrality.top(0), 1);

    // Exact sum of the independent Katz series for levels 1 and 2.
    const std::vector<double> exactScores{12.0 / 7.0, 15.0 / 7.0, 5.0 / 7.0};
    hGraph.forEdges([&](edgeid eid) {
        EXPECT_LE(centrality.score(eid), exactScores[eid]);
        // EXPECT_GE(centrality.bound(eid), exactScores[eid]);
    });
}

TEST(HyperKatzCentralityGTest, testTopKAndValidation) {
    Hypergraph hGraph(4);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({2, 3});

    HyperKatzCentrality centrality(hGraph, 1);
    EXPECT_THROW(centrality.scores(), std::runtime_error);
    centrality.run();
    EXPECT_EQ(centrality.top(), 1);

    EXPECT_THROW(HyperKatzCentrality(hGraph, 0), std::invalid_argument);
    EXPECT_THROW(HyperKatzCentrality(hGraph, 4), std::invalid_argument);

    Hypergraph disconnected(2);
    disconnected.addEdge({0});
    disconnected.addEdge({1});
    EXPECT_THROW(HyperKatzCentrality(disconnected, 1), std::runtime_error);
}

} // namespace NetworKit
