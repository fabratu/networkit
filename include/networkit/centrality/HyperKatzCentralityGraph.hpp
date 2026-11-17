#ifndef NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_GRAPH_HPP_
#define NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_GRAPH_HPP_

#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/SLevelGraphContainer.hpp>

namespace NetworKit {

/**
 * Computes top-k hyperedge Katz scores by running one Katz series per nonempty s-level Graph.
 * Each level uses alpha_s = 1 / (maximumDegree_s + 1), as in DynKatzCentrality.
 */
class HyperKatzCentralityGraph final : public Algorithm {
public:
    HyperKatzCentralityGraph(const Hypergraph &hGraph, count k, bool groupOnly = false,
                             double tolerance = 1e-9);

    void run() override;

    const Vector &scores() const;
    double score(edgeid eid) const;
    std::vector<std::pair<edgeid, double>> ranking() const;
    edgeid top(count n = 0) const;
    double bound(edgeid eid) const;
    bool areDistinguished(edgeid eid1, edgeid eid2) const;
    double getAlpha(count s) const;

    count iterationReached{0};

private:
    void doIteration();
    bool checkConvergence();
    bool areSufficientlyRanked(edgeid high, edgeid low) const;

    const Hypergraph &hGraph;
    SLevelGraphContainer graphs;
    const count k;
    const bool groupOnly;
    const double rankTolerance;

    std::vector<count> levelNumbers;
    std::vector<double> levelAlphas;
    std::vector<double> levelUpperFactors;
    std::vector<double> alphaByLevel;
    std::vector<std::vector<double>> currentPaths;
    std::vector<std::vector<double>> nextPaths;
    std::vector<double> alphaPowers;
    Vector msLowerBound;
    Vector msUpperBound;
    std::vector<edgeid> activeRanking;
};

} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_HYPER_KATZ_CENTRALITY_GRAPH_HPP_
