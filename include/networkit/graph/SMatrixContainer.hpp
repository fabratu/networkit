#ifndef NETWORKIT_GRAPH_S_MATRIX_CONTAINER_HPP_
#define NETWORKIT_GRAPH_S_MATRIX_CONTAINER_HPP_

#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/** Selects which family of s-matrices a container builds. */
enum class SMatrixType { Level, Delta };

/**
 * Stores s-level adjacency matrices or s-delta matrices of a hypergraph.
 *
 * Levels are built from one through one level beyond the maximum hyperedge intersection. Equal
 * s-level matrices are stored only once. Each non-empty s-delta matrix is stored separately, while
 * delta levels without interactions point to one shared empty matrix. The referenced hypergraph
 * must outlive this container.
 */
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
     * @param type Whether to build s-level or s-delta matrices.
     */
    void build(SMatrixType type = SMatrixType::Level);

    /** Removes all matrices and level mappings from the container. */
    void reset() noexcept;

    /**
     * Returns the adjacency matrix for level @a s.
     *
     * @throws std::invalid_argument If @a s is zero.
     * @throws std::out_of_range If @a s is greater than getMaxLevel().
     * @throws std::runtime_error If build() has not been called since construction or reset().
     */
    const CSRMatrix &getMatrix(count s) const;

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
            const CSRMatrix &matrix = matrices[levelToMatrix[s - 1]];
            if (matrix.nnz() != 0)
                handle(s, matrix);
        }
    }

private:
    const Hypergraph &hGraph;
    std::vector<CSRMatrix> matrices;
    std::vector<index> levelToMatrix;
    SMatrixType type{SMatrixType::Level};
    bool built{false};
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_S_MATRIX_CONTAINER_HPP_
