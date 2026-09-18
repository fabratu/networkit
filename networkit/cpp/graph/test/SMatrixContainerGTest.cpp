#include <gtest/gtest.h>

#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

TEST(SMatrixContainerGTest, testBuildAndLevelLookup) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    SMatrixContainer container(hGraph);
    EXPECT_FALSE(container.isBuilt());
    EXPECT_EQ(container.numberOfMatrices(), 0);

    container.build();

    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getType(), SMatrixType::Level);
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

TEST(SMatrixContainerGTest, testDeduplicatesLevels) {
    Hypergraph hGraph(2);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1});

    SMatrixContainer container(hGraph);
    container.build();

    EXPECT_EQ(container.getMaxLevel(), 3);
    EXPECT_EQ(container.numberOfMatrices(), 2);
    EXPECT_EQ(container.getMatrix(1).nnz(), 2);
    EXPECT_EQ(container.getMatrix(3).nnz(), 0);
    EXPECT_EQ(&container.getMatrix(1), &container.getMatrix(2));
}

TEST(SMatrixContainerGTest, testEmptyHypergraph) {
    Hypergraph hGraph;
    SMatrixContainer container(hGraph);
    container.build();

    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getMaxLevel(), 1);
    EXPECT_EQ(container.numberOfMatrices(), 1);
    EXPECT_EQ(container.getMatrix(1).numberOfRows(), 0);
    EXPECT_EQ(container.getMatrix(1).nnz(), 0);
    EXPECT_THROW(container.getMatrix(2), std::out_of_range);
}

TEST(SMatrixContainerGTest, testBuildDeltaMatrices) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1});
    hGraph.addEdge({2, 3, 4, 5});
    hGraph.addEdge({2, 3, 4, 5});

    SMatrixContainer container(hGraph);
    container.build(SMatrixType::Delta);

    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getType(), SMatrixType::Delta);
    EXPECT_EQ(container.getMaxLevel(), 5);
    EXPECT_EQ(container.numberOfMatrices(), 3);

    EXPECT_EQ(container.getMatrix(1).nnz(), 0);
    EXPECT_EQ(container.getMatrix(2).nnz(), 2);
    EXPECT_DOUBLE_EQ(container.getMatrix(2)(0, 1), 1.0);
    EXPECT_EQ(container.getMatrix(3).nnz(), 0);
    EXPECT_EQ(container.getMatrix(4).nnz(), 2);
    EXPECT_DOUBLE_EQ(container.getMatrix(4)(2, 3), 1.0);
    EXPECT_EQ(container.getMatrix(5).nnz(), 0);

    EXPECT_EQ(&container.getMatrix(1), &container.getMatrix(3));
    EXPECT_EQ(&container.getMatrix(3), &container.getMatrix(5));
    EXPECT_NE(&container.getMatrix(2), &container.getMatrix(4));
}

TEST(SMatrixContainerGTest, testForLevels) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1});
    hGraph.addEdge({2, 3, 4, 5});
    hGraph.addEdge({2, 3, 4, 5});

    SMatrixContainer container(hGraph);
    EXPECT_THROW(container.forLevels([](count, const CSRMatrix &) {}), std::runtime_error);

    container.build(SMatrixType::Level);
    std::vector<count> levelLevels;
    std::vector<const CSRMatrix *> levelMatrices;
    container.forLevels([&](count s, const CSRMatrix &matrix) {
        levelLevels.push_back(s);
        levelMatrices.push_back(&matrix);
    });

    EXPECT_EQ(levelLevels, (std::vector<count>{1, 2, 3, 4}));
    ASSERT_EQ(levelMatrices.size(), 4);
    EXPECT_EQ(levelMatrices[0], levelMatrices[1]);
    EXPECT_NE(levelMatrices[1], levelMatrices[2]);
    EXPECT_EQ(levelMatrices[2], levelMatrices[3]);

    container.build(SMatrixType::Delta);
    std::vector<count> deltaLevels;
    container.forLevels([&](count s, const CSRMatrix &matrix) { deltaLevels.push_back(s); });

    EXPECT_EQ(deltaLevels, (std::vector<count>{2, 4}));
}

TEST(SMatrixContainerGTest, testResetAndRebuild) {
    Hypergraph hGraph(2);
    hGraph.addEdge({0});

    SMatrixContainer container(hGraph);
    container.build();
    ASSERT_TRUE(container.isBuilt());

    container.reset();
    EXPECT_FALSE(container.isBuilt());
    EXPECT_EQ(container.numberOfMatrices(), 0);
    EXPECT_EQ(container.getMaxLevel(), 0);
    EXPECT_THROW(container.getMatrix(1), std::runtime_error);
    EXPECT_THROW(container.getType(), std::runtime_error);

    hGraph.addEdge({0, 1});
    container.build();
    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getMatrix(1).numberOfRows(), 2);
    EXPECT_DOUBLE_EQ(container.getMatrix(1)(0, 1), 1.0);
}

} // namespace NetworKit
