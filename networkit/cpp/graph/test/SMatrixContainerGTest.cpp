#include <gtest/gtest.h>

#include <networkit/algebraic/ACSRMatrix.hpp>
#include <networkit/algebraic/DCSRMatrix.hpp>
#include <networkit/algebraic/VSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

TEST(SMatrixContainerGTest, testBuildAndLevelLookup) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    SMatrixContainer<> container(hGraph);
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

    SMatrixContainer<> container(hGraph);
    container.build();

    EXPECT_EQ(container.getMaxLevel(), 3);
    EXPECT_EQ(container.numberOfMatrices(), 2);
    EXPECT_EQ(container.getMatrix(1).nnz(), 2);
    EXPECT_EQ(container.getMatrix(3).nnz(), 0);
    EXPECT_EQ(&container.getMatrix(1), &container.getMatrix(2));
}

TEST(SMatrixContainerGTest, testEmptyHypergraph) {
    Hypergraph hGraph;
    SMatrixContainer<> container(hGraph);
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

    SMatrixContainer<> container(hGraph);
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

    SMatrixContainer<> container(hGraph);
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

    SMatrixContainer<> container(hGraph);
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

TEST(SMatrixContainerGTest, testDCSRStorage) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    SMatrixContainer<DCSRMatrix> container(hGraph);
    container.build();

    const Vector ones(hGraph.upperEdgeIdBound(), 1.0);
    EXPECT_EQ(container.getMatrix(1) * ones, Vector({2.0, 2.0, 2.0, 0.0}));
    EXPECT_EQ(container.getMatrix(2) * ones, Vector({1.0, 1.0, 0.0, 0.0}));
    EXPECT_EQ(container.getMatrix(3) * ones, Vector(4));

    count visitedLevels = 0;
    container.forLevels([&](count, const DCSRMatrix &matrix) {
        EXPECT_EQ(matrix.numberOfRows(), hGraph.upperEdgeIdBound());
        ++visitedLevels;
    });
    EXPECT_EQ(visitedLevels, container.getMaxLevel());

    container.build(SMatrixType::Delta);
    EXPECT_EQ(container.getMatrix(1) * ones, Vector({1.0, 1.0, 2.0, 0.0}));
    EXPECT_EQ(container.getMatrix(2) * ones, Vector({1.0, 1.0, 0.0, 0.0}));
    EXPECT_EQ(container.getMatrix(3) * ones, Vector(4));
}

TEST(SMatrixContainerGTest, testPatternOnlyDCSRDeltaStorage) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    SMatrixContainer<DCSRMatrix> dcsr(hGraph);
    SMatrixContainer<ACSRMatrix> acsr(hGraph);
    std::vector<double> dcsrAlphas;
    std::vector<double> acsrAlphas;
    dcsr.build(SMatrixType::Delta, nullptr, &dcsrAlphas, true);
    acsr.build(SMatrixType::Delta, nullptr, &acsrAlphas, true);

    EXPECT_EQ(dcsrAlphas, acsrAlphas);
    EXPECT_EQ(dcsr.getMaxLevel(), acsr.getMaxLevel());
    const Vector input({1.0, 2.0, 3.0, 4.0});
    for (count level = 1; level <= dcsr.getMaxLevel(); ++level) {
        const auto &dcsrMatrix = dcsr.getMatrix(level);
        const auto &acsrMatrix = acsr.getMatrix(level);
        EXPECT_TRUE(dcsrMatrix.hasUnitValues());
        EXPECT_EQ(dcsrMatrix.nnz(), acsrMatrix.nnz());
        Vector acsrResult(hGraph.upperEdgeIdBound(), 0.0);
        acsrMatrix.addUnitValuesInto(input, acsrResult);
        EXPECT_EQ(dcsrMatrix * input, acsrResult);
        Vector accumulated(hGraph.upperEdgeIdBound(), 3.0);
        dcsrMatrix.addUnitValuesInto(input, accumulated);
        EXPECT_EQ(accumulated, acsrResult + Vector(hGraph.upperEdgeIdBound(), 3.0));
    }
}

TEST(SMatrixContainerGTest, testPatternOnlyDCSRLevelStorage) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    SMatrixContainer<DCSRMatrix> pattern(hGraph);
    SMatrixContainer<CSRMatrix> valued(hGraph);
    std::vector<double> patternAlphas;
    std::vector<double> valuedAlphas;
    pattern.build(SMatrixType::Level, nullptr, &patternAlphas, true);
    valued.build(SMatrixType::Level, nullptr, &valuedAlphas);

    EXPECT_EQ(patternAlphas, valuedAlphas);
    EXPECT_EQ(pattern.getMaxLevel(), valued.getMaxLevel());
    const Vector input({1.0, 2.0, 3.0, 4.0});
    for (count level = 1; level < pattern.getMaxLevel(); ++level) {
        const auto &matrix = pattern.getMatrix(level);
        EXPECT_TRUE(matrix.hasUnitValues());
        EXPECT_EQ(matrix.nnz(), valued.getMatrix(level).nnz());
        EXPECT_EQ(matrix * input, valued.getMatrix(level) * input);

        Vector output(hGraph.upperEdgeIdBound(), 7.0);
        matrix.multiplyUnitValuesInto(input, output);
        EXPECT_EQ(output, valued.getMatrix(level) * input);
    }
    EXPECT_EQ(pattern.getMatrix(pattern.getMaxLevel()).nnz(), 0);

    Hypergraph repeated(2);
    repeated.addEdge({0, 1});
    repeated.addEdge({0, 1});
    SMatrixContainer<DCSRMatrix> repeatedPattern(repeated);
    repeatedPattern.build(SMatrixType::Level, nullptr, nullptr, true);
    EXPECT_EQ(&repeatedPattern.getMatrix(1), &repeatedPattern.getMatrix(2));
    EXPECT_EQ(repeatedPattern.getMatrix(3).nnz(), 0);
}

TEST(SMatrixContainerGTest, testACSRDeltaStorage) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 1});
    hGraph.addEdge({2, 3, 4, 5});
    hGraph.addEdge({2, 3, 4, 5});

    SMatrixContainer<ACSRMatrix> container(hGraph);
    EXPECT_THROW(container.build(SMatrixType::Level), std::invalid_argument);
    EXPECT_FALSE(container.isBuilt());

    container.build(SMatrixType::Delta);
    EXPECT_TRUE(container.isBuilt());
    EXPECT_EQ(container.getType(), SMatrixType::Delta);
    EXPECT_EQ(container.getMaxLevel(), 5);
    EXPECT_EQ(container.numberOfMatrices(), 5);

    for (count s = 1; s < container.getMaxLevel(); ++s)
        EXPECT_NE(&container.getMatrix(s), &container.getMatrix(s + 1));

    count visitedLevels = 0;
    container.forLevels([&](count, const ACSRMatrix &) { ++visitedLevels; });
    EXPECT_EQ(visitedLevels, container.getMaxLevel());
}

TEST(SMatrixContainerGTest, testComputesLevelAlphasDuringDeltaBuild) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    std::vector<double> alphas;
    SMatrixContainer<ACSRMatrix> container(hGraph);
    container.build(SMatrixType::Delta, nullptr, &alphas);

    ASSERT_EQ(alphas.size(), container.getMaxLevel());
    EXPECT_DOUBLE_EQ(alphas[0], 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(alphas[1], 0.5);
    EXPECT_DOUBLE_EQ(alphas[2], 1.0);
}

TEST(SMatrixContainerGTest, testParallelPatternOnlyACSRBuildAndMultiply) {
    Hypergraph hGraph(1);
    for (count edge = 0; edge < 10; ++edge)
        hGraph.addEdge({0});

    SMatrixContainer<ACSRMatrix> container(hGraph);
    VSRMatrix input;
    std::vector<double> alphas;
    container.build(SMatrixType::Delta, &input, &alphas, true);

    ASSERT_GE(alphas.size(), 1);
    EXPECT_DOUBLE_EQ(alphas[0], 0.1);

    VSRMatrix output = input;
    output.reset();
    container.multiplyInto(input, output);
    for (index row = 0; row < hGraph.numberOfEdges(); ++row)
        EXPECT_DOUBLE_EQ(output(row, 0), 9.0);

    const auto &pattern = container.getMatrix(1);
    EXPECT_EQ(pattern.nnz(), 90);
    Vector vectorInput(hGraph.numberOfEdges(), 1.0);
    Vector vectorOutput(hGraph.numberOfEdges(), 2.0);
    pattern.addUnitValuesInto(vectorInput, vectorOutput);
    for (index row = 0; row < hGraph.numberOfEdges(); ++row)
        EXPECT_DOUBLE_EQ(vectorOutput[row], 11.0);
}

TEST(SMatrixContainerGTest, testCompressedACSRPatternPreservesProducts) {
    const std::vector<index> offsets{0, 0, 1, 1, 1, 2, 2};
    const std::vector<index> columns{2, 1};
    ACSRMatrix dense(6, 6, 2, std::vector<index>(offsets), std::vector<index>(columns));
    ACSRMatrix compressed(6, 6, 2, std::vector<index>(offsets),
                          std::vector<index>(columns), false, true);

    ASSERT_FALSE(dense.hasCompressedRows());
    ASSERT_TRUE(compressed.hasCompressedRows());
    EXPECT_LT(compressed.rowIndexBytes(), dense.rowIndexBytes());

    const Vector vectorInput({1.0, 5.0, 3.0, 0.0, 0.0, 0.0});
    Vector denseOutput(6, 2.0);
    Vector compressedOutput(6, 2.0);
    dense.addUnitValuesInto(vectorInput, denseOutput);
    compressed.addUnitValuesInto(vectorInput, compressedOutput);
    EXPECT_EQ(compressedOutput, denseOutput);
    EXPECT_EQ(compressedOutput, Vector({2.0, 5.0, 2.0, 2.0, 7.0, 2.0}));

    VSRMatrix input(6, 3, 2, std::vector<count>(6, 2));
    input(2, 0) = 3.0;
    input(2, 1) = 4.0;
    input(1, 0) = 5.0;
    input(1, 1) = 6.0;
    VSRMatrix denseProduct(6, 3, 2, std::vector<count>(6, 2));
    VSRMatrix compressedProduct = denseProduct;
    dense.multiplyInto(input, denseProduct);
    compressed.multiplyInto(input, compressedProduct);
    for (index row = 0; row < 6; ++row) {
        for (index value = 0; value < 2; ++value)
            EXPECT_DOUBLE_EQ(compressedProduct(row, value), denseProduct(row, value));
    }

    ACSRMatrix empty(6, 6, 1, std::vector<index>(7, 0), std::vector<index>{}, false, true);
    ASSERT_TRUE(empty.hasCompressedRows());
    EXPECT_EQ(empty.rowIndexBytes(), sizeof(index));
    EXPECT_THROW(empty.updateOther(denseProduct), std::logic_error);
}

TEST(SMatrixContainerGTest, testCombinedMultiplyWithCompressedLevel) {
    Hypergraph hGraph(12);
    for (node u = 1; u <= 10; ++u) {
        if (u <= 2)
            hGraph.addEdge({0, u, 11});
        else
            hGraph.addEdge({0, u});
    }

    SMatrixContainer<ACSRMatrix> container(hGraph);
    VSRMatrix input;
    container.build(SMatrixType::Delta, &input, nullptr, true);
    ASSERT_FALSE(container.getMatrix(1).hasCompressedRows());
    ASSERT_TRUE(container.getMatrix(2).hasCompressedRows());

    VSRMatrix combined = input;
    combined.reset();
    container.multiplyInto(input, combined);

    VSRMatrix separate = input;
    separate.reset();
    container.forLevels([&](count, const ACSRMatrix &matrix) {
        matrix.multiplyInto(input, separate);
    });
    for (index row = 0; row < hGraph.numberOfEdges(); ++row) {
        EXPECT_DOUBLE_EQ(combined(row, 0), 9.0);
        EXPECT_DOUBLE_EQ(combined(row, 0), separate(row, 0));
    }
    for (index row = 0; row < 2; ++row) {
        EXPECT_DOUBLE_EQ(combined(row, 1), 1.0);
        EXPECT_DOUBLE_EQ(combined(row, 1), separate(row, 1));
    }
}

TEST(SMatrixContainerGTest, testBuildSetsVSRMatrixWithOneRowPerEdge) {
    Hypergraph hGraph(6);
    hGraph.addEdge({0, 1, 2});
    hGraph.addEdge({0, 1, 3});
    hGraph.addEdge({0, 4});
    hGraph.addEdge({5});

    VSRMatrix vsrMatrix;
    SMatrixContainer<> container(hGraph);
    container.build(SMatrixType::Level, &vsrMatrix);

    // The largest intersection is two, hence the VSR matrix has three columns. Its values are
    // initialized to zero and it contains exactly one result row per hyperedge.
    EXPECT_EQ(vsrMatrix * Vector(3, 1.0), Vector(hGraph.numberOfEdges(), 0.0));
}

TEST(SMatrixContainerGTest, testBuildSetsCompactVSRMatrixAfterEdgeRemoval) {
    Hypergraph hGraph(3);
    hGraph.addEdge({0, 1});
    hGraph.addEdge({0, 2});
    hGraph.addEdge({1, 2});
    hGraph.removeEdge(1);

    VSRMatrix vsrMatrix;
    SMatrixContainer<> container(hGraph);
    container.build(SMatrixType::Delta, &vsrMatrix);

    EXPECT_EQ(vsrMatrix * Vector(2, 1.0), Vector(hGraph.numberOfEdges(), 0.0));
}

} // namespace NetworKit
