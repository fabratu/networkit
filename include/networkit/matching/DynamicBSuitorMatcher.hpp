#ifndef NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_

#include <set>
#include <unordered_map>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/matching/BSuitorMatcher.hpp>

namespace NetworKit {
class DynamicBSuitorMatcher final : public BSuitorMatcher {
    enum class Operation { Insert, Remove };

    struct PairHash
    {
        template <class T1, class T2>
        std::size_t operator() (const std::pair<T1, T2> &v) const
        {
            return std::hash<T1>()(v.first) ^ std::hash<T2>()(v.second) << 1;
        }
    };

    bool numberOfAffectedEquals(int num) {
        return std::count_if(affected.begin(), affected.end(), [](bool aff) { return aff; }) == num;
    }

    bool isBetterMatch(node u, node v, edgeweight ew) const noexcept {
        const auto currentMatch = Suitors.at(u)->min;
        bool isBetterMatch = currentMatch.id == none || currentMatch.weight < ew
               || (currentMatch.weight == ew && v < currentMatch.id);
        // INFO("Node ", u, " Min: ", currentMatch.id, " Min weight: ", currentMatch.weight, " ", ew, " is better: ", isBetterMatch);
        return isBetterMatch;
    }

    void processEdgeInsertion(const WeightedEdge &edge);
    void processEdgeRemoval(const Edge &edge);

    // New working code
    void findAffectedNodes(node u, node v, Operation op);
    void updateAffectedNodes(uint8_t batchId);

    // Other approach
    void trackUpdatePath(size_t batchId, node start, bool recursiveCall = false);
    void processEdgeInsertionNew(const WeightedEdge &edge);
    void processEdgeRemovalNew(const Edge &edge);
    // void updateAffectedNodes();

    // // Original code
    // void findAffectedNodes(node u, node v, Operation op);
    // void updateAffectedNodes();

    // // Using S-invariant consequently + update only if better choice
    // void findAffectedNodes2(node u, node v, Operation op);
    // void updateAffectedNodes2();

    // // Using T-invariant + S-invariant consequently
    // void findAffectedNodes3(node u, node v, Operation op);
    // void updateAffectedNodes3();

    // // Using T-invariant + S-invariant consequently + count visits
    // void findAffectedNodes4(node u, node v, Operation op);
    // void updateAffectedNodes4();

    // // Using S-invariant consequently + update only if better choice + fix loose end
    // void findAffectedNodes5(node u, node v, Operation op);
    // void updateAffectedNodes5();

    // // Using S-invariant consequently + update only if better choice + search for unsaturated
    // loose
    // // ends
    // void findAffectedNodes6(node u, node v, Operation op);
    // void updateAffectedNodes6();

    // // Using S-invariant consequently + update only if better choice + search for loose ends +
    // allow
    // // circles for removal
    // void findAffectedNodes7(node u, node v, Operation op);
    // void updateAffectedNodes7(Operation op);

    // // Using S-invariant consequently + update only if better choice + search for saturated loose
    // // ends in paths
    // void findAffectedNodes8(node u, node v, Operation op);
    // void updateAffectedNodes8();

public:
    DynamicBSuitorMatcher(const Graph &G, const std::vector<count> &b) : BSuitorMatcher(G, b) {
        // affected.resize(G.upperNodeIdBound());
        affectedNodes.reserve(G.numberOfNodes());
    }
    DynamicBSuitorMatcher(const Graph &G, count b = 1) : BSuitorMatcher(G, b) {
        // affected.resize(G.upperNodeIdBound());
        affectedNodes.reserve(G.numberOfNodes());
    }
    DynamicBSuitorMatcher(const Graph &G, const std::string &path) : BSuitorMatcher(G, b) {
        // affected.resize(G.upperNodeIdBound());
        affectedNodes.reserve(G.numberOfNodes());
    }

    std::vector<Edge> edgeBatch;
    std::vector<node> affectedNodes;
    std::vector<bool> affected;
    count affectedNodesPerRun;
    // TODO: support bigger batch-sizes
    std::unordered_map<std::pair<int,int>,size_t,PairHash> batchTracker;

    void addEdge(WeightedEdge &edge) {
        std::vector<WeightedEdge> edges{edge};
        addEdges(edges);
    }
    void addEdges(std::vector<WeightedEdge> &edges, bool sort = true);

    void removeEdge(Edge &edge) {
        std::vector<Edge> edges{edge};
        removeEdges(edges);
    }
    void removeEdges(std::vector<Edge> &edges);

    count getNumberOfAffected();
};

} // namespace NetworKit
#endif // NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
