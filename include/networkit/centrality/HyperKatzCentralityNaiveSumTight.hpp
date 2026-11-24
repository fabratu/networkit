#ifndef NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_NAIVE_SUM_TIGHT_HPP_
#define NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_NAIVE_SUM_TIGHT_HPP_

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
 * Each non-empty s-level matrix A_s contributes an independent Katz series with the
 * row-normalized operator D_s A_s, where D_s[u,u] = 1 / (degree_s(u) + 1).
 * The score is the sum of the series, excluding their zero-length terms.
 * The stopping test uses a degree-local tail factor. This factor is not a
 * certified upper bound for arbitrary irregular s-level graphs.
 *
 * @see https://doi.org/10.4230/LIPIcs.ESA.2018.42
 */
class HyperKatzCentralityNaiveSumTight final : public Algorithm {
public:
    /**
     * Constructs the algorithm for @a hGraph. Each hyperedge in a non-empty level receives
     * a damping factor 1 / (degree_s(eid) + 1).
     *
     * @param hGraph Input hypergraph.
     * @param k Number of highest-ranked hyperedges to identify.
     * @param groupOnly Whether only membership in the top-k, rather than its order, is required.
     * @param tolerance Ranking tolerance used by the convergence test.
     */
    HyperKatzCentralityNaiveSumTight(const Hypergraph &hGraph, count k, bool groupOnly = false,
                                     double tolerance = 1e-9);

    void run() override;

    /** Returns scores indexed by hyperedge id. */
    const Vector &scores() const;

    /** Returns the score of hyperedge @a eid. */
    double score(edgeid eid) const;

    /** Returns all existing hyperedges sorted by decreasing score. */
    std::vector<std::pair<edgeid, double>> ranking() const;

    /** Returns the hyperedge at position @a n in the computed top-k ranking. */
    edgeid top(count n = 0) const;

    /** Returns the degree-local convergence estimate of hyperedge @a eid. */
    double bound(edgeid eid) const;

    /** Returns whether the bounds establish an order between two hyperedges. */
    bool areDistinguished(edgeid eid1, edgeid eid2) const;

    /**
     * Returns the damping factor selected for hyperedge @a eid in level @a s.
     *
     * @throws std::invalid_argument If the level's matrix is empty.
     * @throws std::out_of_range If @a s is outside the container's levels.
     */
    double getAlpha(count s, edgeid eid) const;

    /** Number of completed approximation iterations. */
    count iterationReached{0};

private:
    void doIteration();
    bool checkConvergence();
    bool areSufficientlyRanked(edgeid high, edgeid low) const;

    const Hypergraph &hGraph;
    SMatrixContainer<DCSRMatrix> matrices;
    const count k;
    const bool groupOnly;
    const double rankTolerance;

    std::vector<const DCSRMatrix *> levelMatrices;
    std::vector<Vector> levelAlphas;
    std::vector<Vector> levelUpperFactors;
    std::vector<index> levelIndexByNumber;
    std::vector<Vector> currentPaths;
    std::vector<Vector> nextPaths;
    Vector msLowerBound;
    Vector msUpperBound;
    std::vector<edgeid> activeRanking;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_NAIVE_SUM_TIGHT_HPP_
