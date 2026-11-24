#ifndef NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_DELTA_IN_PLACE_TIGHT_HPP_
#define NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_DELTA_IN_PLACE_TIGHT_HPP_

#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/VSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

/**
 * Computes a top-k Katz ranking of the hyperedges of a hypergraph.
 *
 * Uses the delta-matrix traversal of HyperKatzCentralityDeltaInPlace with the
 * per-hyperedge, per-level damping and bound factors of
 * HyperKatzCentralityNaiveSumTight.
 *
 * @see https://doi.org/10.4230/LIPIcs.ESA.2018.42
 */
class HyperKatzCentralityDeltaInPlaceTight final : public Algorithm {
public:
    /**
     * Constructs the algorithm for @a hGraph. Each hyperedge in level s receives
     * damping factor 1 / (degree_s(eid) + 1).
     *
     * @param hGraph Input hypergraph.
     * @param k Number of highest-ranked hyperedges to identify.
     * @param groupOnly Whether only membership in the top-k, rather than its order, is required.
     * @param tolerance Ranking tolerance used by the convergence test.
     */
    HyperKatzCentralityDeltaInPlaceTight(const Hypergraph &hGraph, count k, bool groupOnly = false,
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

    /** Returns the current ranking-bound estimate of hyperedge @a eid. */
    double bound(edgeid eid) const;

    /** Returns whether the bounds establish an order between two hyperedges. */
    bool areDistinguished(edgeid eid1, edgeid eid2) const;

    /** Number of completed approximation iterations. */
    count iterationReached{0};

private:
    void doIteration();
    bool checkConvergence();
    bool areSufficientlyRanked(edgeid high, edgeid low) const;

    const Hypergraph &hGraph;
    SMatrixContainer<ACSRMatrix> matrices;
    const count k;
    const bool groupOnly;
    const double rankTolerance;

    VSRMatrix currentPaths;
    VSRMatrix nextPaths;
    VSRMatrix levelAlphas;
    VSRMatrix levelUpperFactors;
    Vector msLowerBound;
    Vector msUpperBound;
    std::vector<edgeid> activeRanking;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_DELTA_IN_PLACE_TIGHT_HPP_
