#ifndef NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_NAIVE_SUM_PER_LEVEL_EPS_HPP_
#define NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_NAIVE_SUM_PER_LEVEL_EPS_HPP_

#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

/**
 * Computes a top-k Katz ranking of the hyperedges of a hypergraph.
 *
 * The algorithm adapts the bound-based approximation scheme from "Scalable Katz Ranking
 * Computation in Large Static and Dynamic Graphs". Its effective adjacency operator is the sum
 * of the independent Katz series of all non-empty s-level adjacency matrices. Each level s uses
 * its own damping factor alpha_s = 1 / (d_s + 1), where d_s is the maximum degree of its matrix.
 * Each iteration applies every s-level matrix to its own path vector with an SpMV and sums the
 * resulting score and bound contributions.
 *
 * @see https://doi.org/10.4230/LIPIcs.ESA.2018.42
 */
class HyperKatzCentralityNaiveSumPerLevelEps final : public Algorithm {
public:
    /**
     * Constructs the algorithm for @a hGraph. Each non-empty level receives a damping factor
     * 1 / (d_s + 1), where d_s is the maximum row sum of that level's matrix.
     *
     * @param hGraph Input hypergraph.
     * @param k Number of highest-ranked hyperedges to identify.
     * @param groupOnly Whether only membership in the top-k, rather than its order, is required.
     * @param tolerance Ranking tolerance used by the convergence test.
     */
    HyperKatzCentralityNaiveSumPerLevelEps(const Hypergraph &hGraph, count k,
                                           bool groupOnly = false, double tolerance = 1e-9);

    void run() override;

    /** Returns scores indexed by hyperedge id. */
    const Vector &scores() const;

    /** Returns the score of hyperedge @a eid. */
    double score(edgeid eid) const;

    /** Returns all existing hyperedges sorted by decreasing score. */
    std::vector<std::pair<edgeid, double>> ranking() const;

    /** Returns the hyperedge at position @a n in the computed top-k ranking. */
    edgeid top(count n = 0) const;

    // /** Returns the upper Katz bound of hyperedge @a eid. */
    double bound(edgeid eid) const;

    /** Returns whether the bounds establish an order between two hyperedges. */
    bool areDistinguished(edgeid eid1, edgeid eid2) const;

    /** Number of completed approximation iterations. */
    count iterationReached{0};

private:
    void doIteration();
    bool checkConvergence();
    bool checkGlobalConvergence();
    bool areSufficientlyRanked(edgeid high, edgeid low) const;

    const Hypergraph &hGraph;
    SMatrixContainer matrices;
    const count k;
    const bool groupOnly;
    const double rankTolerance;

    std::vector<const CSRMatrix *> levelMatrices;
    std::vector<count> levelMaxDegrees;
    std::vector<double> levelAlphas;
    std::vector<double> levelTolerances;
    std::vector<Vector> currentPaths;
    std::vector<Vector> lowerBound;
    std::vector<Vector> upperBound;
    Vector msLowerBound;
    Vector msUpperBound;
    std::vector<edgeid> activeRanking;
    count activeLevel;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_NAIVE_SUM_PER_LEVEL_EPS_HPP_
