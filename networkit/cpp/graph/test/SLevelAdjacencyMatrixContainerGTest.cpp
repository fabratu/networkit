#include <gtest/gtest.h>

#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/SLevelAdjacencyMatrixContainer.hpp>

namespace NetworKit {

TEST(SLevelAdjacencyMatrixContainerGTest, testBuildAndLevelLookup) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    SLevelAdjacencyMatrixContainer container(hGraph);
    EXPECT_FALSE(container.isBuilt());
    EXPECT_EQ(container.numberOfMatrices(), 0);

    container.build();

    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getMaxLevel(), 3);
    EXPECT_EQ(container.numberOfMatrices(), 3);
    EXPECT_EQ(container.getMatrix(1).nnz(), 6);
    EXPECT_EQ(container.getMatrix(2).nnz(), 2);
    EXPECT_EQ(container.getMatrix(3).nnz(), 0);
    EXPECT_DOUBLE_EQ(container.getMatrix(1)(0, 2), 1.0);
    EXPECT_DOUBLE_EQ(container.getMatrix(2)(0, 2), 0.0);

    EXPECT_THROW(container.getMatrix(0), std::invalid_argument);
    EXPECT_THROW(container.getMatrix(4), std::out_of_range);
}

TEST(SLevelAdjacencyMatrixContainerGTest, testDeduplicatesLevels) {
    Hypergraph hGraph(2);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1});

    SLevelAdjacencyMatrixContainer container(hGraph);
    container.build();

    EXPECT_EQ(container.getMaxLevel(), 3);
    EXPECT_EQ(container.numberOfMatrices(), 2);
    EXPECT_EQ(container.getMatrix(1).nnz(), 2);
    EXPECT_EQ(container.getMatrix(3).nnz(), 0);
    EXPECT_EQ(&container.getMatrix(1), &container.getMatrix(2));
}

TEST(SLevelAdjacencyMatrixContainerGTest, testEmptyHypergraph) {
    Hypergraph hGraph;
    SLevelAdjacencyMatrixContainer container(hGraph);
    container.build();

    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getMaxLevel(), 1);
    EXPECT_EQ(container.numberOfMatrices(), 1);
    EXPECT_EQ(container.getMatrix(1).numberOfRows(), 0);
    EXPECT_EQ(container.getMatrix(1).nnz(), 0);
    EXPECT_THROW(container.getMatrix(2), std::out_of_range);
}

TEST(SLevelAdjacencyMatrixContainerGTest, testResetAndRebuild) {
    Hypergraph hGraph(2);
    hGraph.addEdge({0});

    SLevelAdjacencyMatrixContainer container(hGraph);
    container.build();
    ASSERT_TRUE(container.isBuilt());

    container.reset();
    EXPECT_FALSE(container.isBuilt());
    EXPECT_EQ(container.numberOfMatrices(), 0);
    EXPECT_EQ(container.getMaxLevel(), 0);
    EXPECT_THROW(container.getMatrix(1), std::runtime_error);

    hGraph.addEdge({0, 1});
    container.build();
    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getMatrix(1).numberOfRows(), 2);
    EXPECT_DOUBLE_EQ(container.getMatrix(1)(0, 1), 1.0);
}

} // namespace NetworKit
