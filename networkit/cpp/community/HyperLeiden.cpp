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
            std::vector<Aux::HTCustodian> edgeCommunityMemberships(G->upperEdgeIdBound());
            std::vector<Aux::HTCustodian> edgeCommunityVolumes(G->upperEdgeIdBound());
            // TODO: implement
            initializeMemberships(*G, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);

            // Greedy Move Phase
            greedyMovePhase(*G, communityMemberships, communitySizes, edgeCommunityMemberships,
                            edgeCommunityVolumes);

            // Refine disconnected communities

            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                std::vector<count> tmpCommunityMemberships(G->upperNodeIdBound(), 0);
                std::vector<count> tmpCommunitySizes(G->upperNodeIdBound(), 0);
                std::vector<Aux::HTCustodian> tmpEdgeCommunityMemberships(G->upperEdgeIdBound());
                std::vector<Aux::HTCustodian> tmpEdgeCommunityVolumes(G->upperEdgeIdBound());
                refineDisconnected(*G, tmpCommunityMemberships, tmpCommunityMemberships,
                                   tmpEdgeCommunityMemberships, tmpEdgeCommunityVolumes,
                                   communityMemberships);
            }

            // Aggregate hypergraph
            currentG = aggregateHypergraph(*G, communityMemberships);
        } else {
            // Initialize memberships + sizes
            std::vector<count> communityMemberships(currentG.upperNodeIdBound(), 0);
            std::vector<count> communitySizes(currentG.upperNodeIdBound(), 0);
            std::vector<Aux::HTCustodian> edgeCommunityMemberships(currentG.upperEdgeIdBound());
            std::vector<Aux::HTCustodian> edgeCommunityVolumes(currentG.upperEdgeIdBound());
            initializeMemberships(currentG, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);

            // Greedy Move Phase
            greedyMovePhase(currentG, communityMemberships, communitySizes,
                            edgeCommunityMemberships, edgeCommunityVolumes);

            // Refine disconnected communities
            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                std::vector<count> tmpCommunityMemberships(currentG.upperNodeIdBound(), 0);
                std::vector<count> tmpCommunitySizes(currentG.upperNodeIdBound(), 0);
                std::vector<Aux::HTCustodian> tmpEdgeCommunityMemberships(G->upperEdgeIdBound());
                std::vector<Aux::HTCustodian> tmpEdgeCommunityVolumes(G->upperEdgeIdBound());
                refineDisconnected(currentG, tmpCommunityMemberships, tmpCommunityMemberships,
                                   tmpEdgeCommunityMemberships, tmpEdgeCommunityVolumes,
                                   communityMemberships);
            }

            // Aggregate hypergraph
            currentG = aggregateHypergraph(currentG, communityMemberships);
        }
    }
}

void HyperLeiden::greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                                  std::vector<count> &communitySizes,
                                  std::vector<Aux::HTCustodian> &edgeCommunityMemberships,
                                  std::vector<Aux::HTCustodian> &edgeCommunityVolumes) {

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
                getBestCommunity(graph, u, communityMemberships, communitySizes,
                                 edgeCommunityMemberships, edgeCommunityVolumes);
            if (bestCommunity != communityMemberships[u]) {
                updateMemberships(graph, u, bestCommunity, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);
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
                                     std::vector<count> &communityMemberships,
                                     std::vector<count> &communitySizes,
                                     std::vector<Aux::HTCustodian> &edgeCommunityMemberships,
                                     std::vector<Aux::HTCustodian> &edgeCommunityVolumes,
                                     std::vector<count> &referenceCommunityMemberships) {

    count maxIter = 100;

    std::vector<bool> vaff(graph.numberOfNodes(), true);
    for (count l = 0; l < maxIter; l++) {
        graph.parallelForNodes([&](node u) {
            if (!vaff[u])
                return;
            if (communitySizes[communityMemberships[u]] != 1) // Only perform for isolated nodes
                return;
            vaff[u] = false;
            auto [bestCommunity, gain] = getBestCommunity(
                graph, u, communityMemberships, communitySizes, edgeCommunityMemberships,
                edgeCommunityVolumes, referenceCommunityMemberships);
            if (bestCommunity != communityMemberships[u]) {
                if (updateMemberships<true>(graph, u, bestCommunity, communityMemberships,
                                            communitySizes, edgeCommunityMemberships,
                                            edgeCommunityVolumes))
                    graph.forNeighborsOf(u, [&](node w) { vaff[w] = true; });
            }
        });
    }
};

Hypergraph HyperLeiden::aggregateHypergraph(const Hypergraph &graph,
                                            const std::vector<count> &communityMemberships) {
    // nothing here yet
    // Plan:
    // 1: Renumber communities based on prefix sum
    // 2: Create new hypergraph with number of nodes = number of communities
    // 3: Iterate over edges of the original hypergraph
    // 3a: For each edge, if singleton, add to new hypergraph and track the community id (vector of
    // comm ids) 3b: For each edge, if not singleton, create hash of sorted membership vectors (hash
    // combine) and add to a hash database, which maps set-hashes to edge ids. TODO: maybe there is
    // something more efficient here
    // 4: For each hash in the hash database, create a new edge in the new hypergraph and add the
    // corresponding node ids with their combined weights.
    // 5: Create a map from the old node ids to the new node ids (?)

    return Hypergraph();
}

// Helper Functions

void HyperLeiden::initializeMemberships(const Hypergraph &graph,
                                        std::vector<count> &communityMemberships,
                                        std::vector<count> &communitySizes,
                                        std::vector<Aux::HTCustodian> &edgeCommunityMemberships,
                                        std::vector<Aux::HTCustodian> &edgeCommunityVolumes) const {
    // In the beginning each node is its own community
    std::iota(std::begin(communityMemberships), std::end(communityMemberships), 0);
    std::fill(std::begin(communitySizes), std::end(communitySizes), 1);

    graph.parallelForEdges([&](edgeid eId, edgeweight /*/*/) {
        auto handleMemberships = edgeCommunityMemberships[eId].make_handle();
        auto handleVolumes = edgeCommunityVolumes[eId].make_handle();
        auto nodes = graph.nodesOf(eId);
        for (auto &it : nodes) {
            handleMemberships->insert(communityMemberships[it.first], 1);
            handleVolumes->insert(communityMemberships[it.first], it.second);
        }
    });
}

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
                              const std::vector<count> &communitySizes,
                              const std::vector<Aux::HTCustodian> &edgeCommunityMemberships,
                              const std::vector<Aux::HTCustodian> &edgeCommunityVolumes,
                              const std::vector<count> &referenceCommunityMemberships) const {
    auto [oldCommunity, oldSize] =
        std::make_pair(communityMemberships[v], communitySizes[communityMemberships[v]]);
    auto [bestCommunity, bestGain] = std::make_pair(0, 0.0);
    auto neighboringCommunities =
        gatherNeighboringCommunities(graph, v, communityMemberships, communitySizes);

    for (const auto &[c, size] : neighboringCommunities) {
        double gain = deltaHCPM(graph, v, oldCommunity, c, oldSize, size, communityMemberships,
                                communitySizes, edgeCommunityMemberships, edgeCommunityVolumes);
        if (gain > bestGain) {
            bestGain = gain;
            bestCommunity = c;
        }
    }
    // TODO: maybe default gain not 0.0
    return {bestCommunity, bestGain};
}

double HyperLeiden::deltaHCPM(const Hypergraph &graph, node v, count c1, count c2, count c1Size,
                              count c2Size, const std::vector<count> &communityMemberships,
                              const std::vector<count> &communitySizes,
                              const std::vector<Aux::HTCustodian> &edgeCommunityMemberships,
                              const std::vector<Aux::HTCustodian> &edgeCommunityVolumes) const {
    double delta_n = 0.0;
    double delta_e = 0.0;
    count rho_total = graph.edgeVolume();

    if (c1 != c2) {
        double delta_n = gamma * rho_total * ((2 << c2Size) - (2 << (c1Size - 1)));
        graph.forEdges([&](edgeid eId, edgeweight eWeight) {
            if (graph.hasNode(v, eId)) {
                nodeweight nWeight = graph.getNodeWeightOf(v, eId);
                count jFactor = edgeCommunityVolumes[eId].current_table()->find(c2) + 2 * nWeight;
                count iFactor = edgeCommunityVolumes[eId].current_table()->find(c1) + nWeight;
                count c2eSize = edgeCommunityMemberships[eId].current_table()->find(c2);
                count c1eSize = edgeCommunityMemberships[eId].current_table()->find(c1);
                delta_e += eWeight * ((2 << c2eSize) * jFactor - (2 << (c1eSize - 1)) * iFactor);
            }
        });
    }
    return delta_e - delta_n;
};

std::vector<bool> HyperLeiden::communityExists(const Hypergraph &graph, node v, count bestCommunity,
                                               const std::vector<count> &communityMemberships,
                                               const std::vector<count> &communitySizes) const {
    std::vector<bool> communityExists(graph.upperNodeIdBound(), false);
#pragma omp parallel for
    for (omp_index i = 0; i < communityMemberships.size(); i++) {
        communityExists[communityMemberships[i]] = true;
    }
    return communityExists;
};

void HyperLeiden::renumberCommunities(const Hypergraph &graph,
                                      const std::vector<count> &communitySizes){

};

void HyperLeiden::prefixSum(std::vector<bool> &communityExists) {
    for (auto entry : communityExists) {

        if (entry) {
            entry = true;
        } else {
            entry = false;
        }
    }
};

} // namespace NetworKit
