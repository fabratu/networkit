/*
 * HyperLeiden.cpp
 *
 * Created on: 09.04.2025
 *     Author: Fabian Brandt-Tumescheit
 */

#include <networkit/community/HyperLeiden.hpp>

namespace NetworKit {
HyperLeiden::HyperLeiden(const Hypergraph &hGraph, int numberOfIterations, double gamma,
                         RefinementStrategy refinementStrategy,
                         CoarseningStrategy coarseningStrategy)
    : CommunityDetectionAlgorithm(hGraph), numberOfIterations(numberOfIterations), gamma(gamma),
      refinementStrategy(refinementStrategy), coarseningStrategy(coarseningStrategy){};

void HyperLeiden::run() {

    int maxPasses = 0;
    bool isFirst = true;
    Hypergraph currentG;

    for (int pass = 0; pass < maxPasses; pass++) {

        if (isFirst) {
            // Initialize memberships + sizes
            std::vector<count> communityMemberships(G->upperNodeIdBound(), 0);
            std::vector<count> communitySizes(G->upperNodeIdBound(), 0);
            Aux::HTCustodian edgeCommunityMemberships(G->upperEdgeIdBound());

            // Greedy Move Phase
            greedyMovePhase(*G, communityMemberships, communitySizes, edgeCommunityMemberships);

            // Refine disconnected communities
            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                refineDisconnected(*G, communityMemberships);
            }

            // Aggregate hypergraph
            currentG = aggregateHypergraph(*G, communityMemberships);
        } else {
            // Initialize memberships + sizes
            std::vector<count> communityMemberships(currentG.upperNodeIdBound(), 0);
            std::vector<count> communitySizes(currentG.upperNodeIdBound(), 0);
            Aux::HTCustodian edgeCommunityMemberships(currentG.upperEdgeIdBound());

            // Greedy Move Phase
            greedyMovePhase(currentG, communityMemberships, communitySizes,
                            edgeCommunityMemberships);

            // Refine disconnected communities
            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                refineDisconnected(currentG, communityMemberships);
            }

            // Aggregate hypergraph
            currentG = aggregateHypergraph(currentG, communityMemberships);
        }
    }
}

void HyperLeiden::greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                                  std::vector<count> &communitySizes,
                                  Aux::HTCustodian &edgeCommunityMemberships) {

    count maxIter = 100;

    std::vector<bool> vaff(graph.numberOfNodes(), true);
    double gainPerRound = 0.0;

    for (count l = 0; l < maxIter; l++) {
        double gainPerRound = 0.0;
        graph.parallelForNodes([&](node u) {
            if (!vaff[u])
                return;
            vaff[u] = false;
            auto [bestCommunity, gain] =
                getBestCommunity(graph, u, communityMemberships, communitySizes);
            if (bestCommunity != communityMemberships[u]) {
                updateMemberships(u, bestCommunity, communityMemberships, communitySizes);
                graph.forNeighborsOf(u, [&](node w) { vaff[w] = true; });
#pragma omp atomic
                gainPerRound += gain;
            }
        });

        if (gainPerRound < tolerance) {
            break;
        }
    }
};

void HyperLeiden::refineDisconnected(const Hypergraph &graph,
                                     std::vector<count> &communityMemberships){
    // nothing here yet
};

Hypergraph HyperLeiden::aggregateHypergraph(const Hypergraph &graph,
                                            const std::vector<count> &communityMemberships) {
    // nothing here yet
    return Hypergraph();
}

// Helper Functions

std::unordered_map<count, count>
HyperLeiden::gatherNeighboringCommunities(const Hypergraph &graph, node v,
                                          const std::vector<count> &communityMemberships,
                                          const std::vector<count> &communitySizes) const {
    std::unordered_map<count, count> communities;
    graph.forNeighborsOf(v, [&](node w) {
        if (communityMemberships[v] != communityMemberships[w]) {
            communities.insert({communityMemberships[w], communitySizes[communityMemberships[w]]});
        }
    });
    return communities;
};

std::pair<count, double>
HyperLeiden::getBestCommunity(const Hypergraph &graph, node v,
                              const std::vector<count> &communityMemberships,
                              const std::vector<count> &communitySizes) {
    auto [oldCommunity, oldSize] =
        std::make_pair(communityMemberships[v], communitySizes[communityMemberships[v]]);
    auto [bestCommunity, bestGain] = std::make_pair(0, 0.0);
    auto neighboringCommunities =
        gatherNeighboringCommunities(graph, v, communityMemberships, communitySizes);

    for (const auto &[c, size] : neighboringCommunities) {
        double gain =
            deltaHCPM(oldCommunity, c, oldSize, size, communityMemberships, communitySizes);
        if (gain > bestGain) {
            bestGain = gain;
            bestCommunity = c;
        }
    }
    // TODO: maybe default gain not 0.0
    return {bestCommunity, bestGain};
};

double HyperLeiden::deltaHCPM(count c1, count c2, count c1Size, count c2Size,
                              const std::vector<count> &communityMemberships,
                              const std::vector<count> &communitySizes) const {
    double delta = 0.0;

    if (c1 != c2) {
        delta -= 1.0;
        // delta = (communitySizes[c1] - communitySizes[c2]) / graph.numberOfEdges();
        // delta += (edgeCommunityMemberships[c1] - edgeCommunityMemberships[c2])
        //          / graph.numberOfEdges();
    }
    return delta;
};

} // namespace NetworKit
