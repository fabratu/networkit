#ifndef NETWORKIT_GRAPH_S_LEVEL_ADJACENCY_MATRIX_CONTAINER_HPP_
#define NETWORKIT_GRAPH_S_LEVEL_ADJACENCY_MATRIX_CONTAINER_HPP_

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/**
 * Stores the distinct s-level adjacency matrices of a hypergraph.
 *
 * Levels are built from one through the first level with an empty adjacency matrix. Equal matrices
 * are stored only once. The referenced hypergraph must outlive this container.
 */
class SLevelAdjacencyMatrixContainer final {
public:
    /**
     * Creates an empty container associated with @a hGraph. Call build() to populate it.
     *
     * @param hGraph Hypergraph whose s-level adjacency matrices shall be stored.
     */
    explicit SLevelAdjacencyMatrixContainer(const Hypergraph &hGraph) : hGraph{hGraph} {}

    /** Builds all distinct s-level adjacency matrices for the current hypergraph. */
    void build();

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

private:
    const Hypergraph &hGraph;
    std::vector<CSRMatrix> matrices;
    std::vector<index> levelToMatrix;
    bool built{false};
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_S_LEVEL_ADJACENCY_MATRIX_CONTAINER_HPP_
