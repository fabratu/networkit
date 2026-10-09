#ifndef NETWORKIT_GRAPH_S_MATRIX_CONTAINER_HPP_
#define NETWORKIT_GRAPH_S_MATRIX_CONTAINER_HPP_

#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/ACSRMatrix.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DCSRMatrix.hpp>
#include <networkit/algebraic/VSRMatrix.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/** Selects which family of s-matrices a container builds. */
enum class SMatrixType { Level, Delta };

/**
 * Stores s-level adjacency matrices or s-delta matrices of a hypergraph.
 *
 * Levels are built from one through one level beyond the maximum hyperedge intersection. Equal
 * s-level matrices are stored only once. Each non-empty s-delta matrix is stored separately, while
 * delta levels without interactions point to one shared empty matrix. ACSR storage is supported
 * for delta matrices only; each level is stored separately because its value count equals the
 * level. The referenced hypergraph must outlive this container.
 *
 * @tparam Matrix Matrix representation used for internal storage. CSRMatrix, DCSRMatrix, and
 * ACSRMatrix are supported.
 */
template <typename Matrix = CSRMatrix>
class SMatrixContainer final {
public:
    /**
     * Creates an empty container associated with @a hGraph. Call build() to populate it.
     *
     * @param hGraph Hypergraph whose s-matrices shall be stored.
     */
    explicit SMatrixContainer(const Hypergraph &hGraph) : hGraph{hGraph} {}

    /**
     * Builds the selected family of s-matrices for the current hypergraph.
     *
     * When @a vsrMatrix is provided, it is replaced with a zero-initialized VSR matrix containing
     * one row per hyperedge. The prefix length of a row is the maximum interaction level of the
     * corresponding edge, or one if the edge has no interactions.
     *
     * @param type Whether to build s-level or s-delta matrices.
     * @param vsrMatrix Optional output for the per-edge variable sparse row matrix.
     * @param alphaVector Optional output for one damping factor per level. Each factor is
     *                    <code>1 / (maximumLevelDegree + 1)</code>.
     */
    void build(SMatrixType type = SMatrixType::Level, VSRMatrix *vsrMatrix = nullptr,
               std::vector<double> *alphaVector = nullptr);

    /** Removes all matrices and level mappings from the container. */
    void reset() noexcept;

    /**
     * Returns the adjacency matrix for level @a s.
     *
     * @throws std::invalid_argument If @a s is zero.
     * @throws std::out_of_range If @a s is greater than getMaxLevel().
     * @throws std::runtime_error If build() has not been called since construction or reset().
     */
    const Matrix &getMatrix(count s) const;

    /** Returns the number of matrices physically stored by the container. */
    count numberOfMatrices() const noexcept { return matrices.size(); }

    /** Returns the highest built level, or zero if the container has not been built. */
    count getMaxLevel() const noexcept { return levelToMatrix.size(); }

    /** Returns true if build() has populated the container. */
    bool isBuilt() const noexcept { return built; }

    /**
     * Returns the kind of matrices created by the last build.
     *
     * @throws std::runtime_error If the container has not been built.
     */
    SMatrixType getType() const;

    /**
     * Iterates over all non-empty levels in ascending order. The handler is called with the level
     * and its corresponding matrix as <code>handle(s, matrix)</code>.
     *
     * @param handle Function called for every non-empty level.
     * @throws std::runtime_error If the container has not been built.
     */
    template <typename L>
    void forLevels(L handle) const {
        if (!built)
            throw std::runtime_error("SMatrixContainer has not been built");

        for (count s = 1; s <= getMaxLevel(); ++s) {
            const Matrix &matrix = matrices[levelToMatrix[s - 1]];
            // if (matrix.nnz() != 0)
            handle(s, matrix);
        }
    }

    template <typename L>
    void forLevelsMutable(L handle) {
        if (!built)
            throw std::runtime_error("SMatrixContainer has not been built");

        for (count s = 1; s <= getMaxLevel(); ++s) {
            Matrix &matrix = matrices[levelToMatrix[s - 1]];
            // if (matrix.nnz() != 0)
            handle(s, matrix);
        }
    }

    template <typename L>
    void forLevelsMutableInParrallel(L handle) {
#pragma omp parallel for schedule(guided)
        for (omp_index i = 0; i < static_cast<omp_index>(getMaxLevel()); ++i) {
            // for (count s = 1; s <= getMaxLevel(); ++s) {
            Matrix &matrix = matrices[levelToMatrix[i]];
            // if (matrix.nnz() != 0)
            handle(i, matrix);
        }
    }

    /**
     * Iterates over all non-empty levels in descending order. The handler is called with the level
     * and its corresponding matrix as <code>handle(s, matrix)</code>.
     *
     * @param handle Function called for every non-empty level.
     * @throws std::runtime_error If the container has not been built.
     */
    template <typename L>
    void forLevelsReverse(L handle) const {
        if (!built)
            throw std::runtime_error("SMatrixContainer has not been built");

        for (count s = getMaxLevel(); s >= 1; --s) {
            const Matrix &matrix = matrices[levelToMatrix[s - 1]];
            // if (matrix.nnz() != 0)
            handle(s, matrix);
        }
    }

    /**
     * Iterates over all non-empty levels in descending order until a certain level. The
     * handler is called with the level and its corresponding matrix as <code>handle(s,
     * matrix)</code>.
     *
     * @param handle Function called for every non-empty level.
     * @throws std::runtime_error If the container has not been built.
     */
    template <typename L>
    void forLevelsReverseUntil(index bound, L handle) const {
        if (!built)
            throw std::runtime_error("SMatrixContainer has not been built");

        for (count s = getMaxLevel(); s >= bound; --s) {
            const Matrix &matrix = matrices[levelToMatrix[s - 1]];
            // if (matrix.nnz() != 0)
            handle(s, matrix);
        }
    }

private:
    const Hypergraph &hGraph;
    std::vector<Matrix> matrices;
    std::vector<index> levelToMatrix;
    SMatrixType type{SMatrixType::Level};
    bool built{false};
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_S_MATRIX_CONTAINER_HPP_
