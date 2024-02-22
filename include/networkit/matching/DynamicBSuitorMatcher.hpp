#ifndef NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_

#include <set>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/matching/BSuitorMatcher.hpp>

namespace NetworKit {
class DynamicBSuitorMatcher final : public BSuitorMatcher {
    enum class Operation { Insert, Remove };

    bool numberOfAffectedEquals(int num) {
        return std::count_if(affected.begin(), affected.end(), [](bool aff) { return aff; }) == num;
    }

    bool isBetterMatch(node u, node v, edgeweight ew) const noexcept {
        const auto currentMatch = Suitors.at(u)->min;
        return currentMatch.weight < ew || (currentMatch.weight == ew && v < currentMatch.id);
    }

    void processEdgeInsertion(const WeightedEdge &edge);
    void processEdgeRemoval(const Edge &edge);

    void findAffectedNodes(node u, node v, Operation op, bool skipFirst = false);
    void updateAffectedNodes();

public:
    DynamicBSuitorMatcher(const Graph &G, const std::vector<count> &b) : BSuitorMatcher(G, b) {
        affected.resize(G.upperNodeIdBound());
        affectedNodes.reserve(G.numberOfNodes());
    }
    DynamicBSuitorMatcher(const Graph &G, count b = 1) : BSuitorMatcher(G, b) {
        affected.resize(G.upperNodeIdBound());
        affectedNodes.reserve(G.numberOfNodes());
    }
    DynamicBSuitorMatcher(const Graph &G, const std::string &path) : BSuitorMatcher(G, b) {
        affected.resize(G.upperNodeIdBound());
        affectedNodes.reserve(G.numberOfNodes());
    }

    std::vector<Edge> edgeBatch;
    std::vector<node> affectedNodes;
    std::vector<bool> affected;
    std::set<node> affectedNodesPerRun;

    void addEdge(WeightedEdge &edge) {
        std::vector<WeightedEdge> edges{edge};
        addEdges(edges);
    }
    void addEdges(std::vector<WeightedEdge> &edges);

    void removeEdge(Edge &edge) {
        std::vector<Edge> edges{edge};
        removeEdges(edges);
    }
    void removeEdges(std::vector<Edge> &edges);

    count getNumberOfAffected();
};

} // namespace NetworKit
#endif // NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
