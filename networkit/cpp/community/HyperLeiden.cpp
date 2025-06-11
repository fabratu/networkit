/*
 * HyperLeiden.cpp
 *
 * Created on: 09.04.2025
 *     Author: Fabian Brandt-Tumescheit
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/HyperLeiden.hpp>

namespace NetworKit {
HyperLeiden::HyperLeiden(const Hypergraph &hGraph, int numberOfIterations, double gamma,
                         double tolerance, RefinementStrategy refinementStrategy,
                         CoarseningStrategy coarseningStrategy)
    : CommunityDetectionAlgorithm(hGraph), numberOfIterations(numberOfIterations), gamma(gamma),
      tolerance(tolerance), refinementStrategy(refinementStrategy),
      coarseningStrategy(coarseningStrategy){};

void HyperLeiden::run() {

    bool isFirst = true;
    Hypergraph currentG;
    mappings.clear();

    for (int pass = 0; pass < numberOfIterations; pass++) {

        if (isFirst) {
            // Initialize memberships + sizes
            std::vector<count> communityMemberships(G->upperNodeIdBound(), 0);
            std::vector<count> communitySizes(G->upperNodeIdBound(), 0);
            std::vector<Aux::ParallelHashMap> edgeCommunityMemberships(G->upperEdgeIdBound());
            std::vector<Aux::ParallelHashMap> edgeCommunityVolumes(G->upperEdgeIdBound());
            auto start = std::chrono::high_resolution_clock::now();
            initializeMemberships(*G, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);
            auto end = std::chrono::high_resolution_clock::now();
            double duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("Initialization time: ", duration, " seconds");
            // Greedy Move Phase
            start = std::chrono::high_resolution_clock::now();
            greedyMovePhase(*G, communityMemberships, communitySizes, edgeCommunityMemberships,
                            edgeCommunityVolumes);
            end = std::chrono::high_resolution_clock::now();
            duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("Greedy move phase time: ", duration, " seconds");

            // Refine disconnected communities
            start = std::chrono::high_resolution_clock::now();
            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                std::vector<count> tmpCommunityMemberships(G->upperNodeIdBound(), 0);
                std::vector<count> tmpCommunitySizes(G->upperNodeIdBound(), 0);
                std::vector<Aux::ParallelHashMap> tmpEdgeCommunityMemberships(
                    G->upperEdgeIdBound());
                std::vector<Aux::ParallelHashMap> tmpEdgeCommunityVolumes(G->upperEdgeIdBound());
                refineDisconnected(*G, tmpCommunityMemberships, tmpCommunityMemberships,
                                   tmpEdgeCommunityMemberships, tmpEdgeCommunityVolumes,
                                   communityMemberships);
            }
            end = std::chrono::high_resolution_clock::now();
            duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("Refinement time: ", duration, " seconds");
            // Aggregate hypergraph
            start = std::chrono::high_resolution_clock::now();
            currentG = aggregateHypergraph(*G, communityMemberships, communitySizes);
            end = std::chrono::high_resolution_clock::now();
            duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("Aggregation time: ", duration, " seconds");
            // Create mapping
            mappings.push_back(createMapping(communityMemberships, communitySizes));

            isFirst = false;
        } else {
            // Initialize memberships + sizes
            std::vector<count> communityMemberships(currentG.upperNodeIdBound(), 0);
            std::vector<count> communitySizes(currentG.upperNodeIdBound(), 0);
            std::vector<Aux::ParallelHashMap> edgeCommunityMemberships(currentG.upperEdgeIdBound());
            std::vector<Aux::ParallelHashMap> edgeCommunityVolumes(currentG.upperEdgeIdBound());
            initializeMemberships(currentG, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);

            // Greedy Move Phase
            greedyMovePhase(currentG, communityMemberships, communitySizes,
                            edgeCommunityMemberships, edgeCommunityVolumes);

            // Refine disconnected communities
            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                std::vector<count> tmpCommunityMemberships(currentG.upperNodeIdBound(), 0);
                std::vector<count> tmpCommunitySizes(currentG.upperNodeIdBound(), 0);
                std::vector<Aux::ParallelHashMap> tmpEdgeCommunityMemberships(
                    G->upperEdgeIdBound());
                std::vector<Aux::ParallelHashMap> tmpEdgeCommunityVolumes(G->upperEdgeIdBound());
                refineDisconnected(currentG, tmpCommunityMemberships, tmpCommunityMemberships,
                                   tmpEdgeCommunityMemberships, tmpEdgeCommunityVolumes,
                                   communityMemberships);
            }

            // Aggregate hypergraph
            currentG = aggregateHypergraph(currentG, communityMemberships, communitySizes);

            // Create mapping
            mappings.push_back(createMapping(communityMemberships, communitySizes));
        }
    }

    // TODO: unroll mappings to get final community memberships
    flattenPartition();

    hasRun = true;
}

void HyperLeiden::greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                                  std::vector<count> &communitySizes,
                                  std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                                  std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) {

    count maxIter = 100;
    count rho_total = graph.edgeVolume();

    std::vector<bool> vaff(graph.numberOfNodes(), true);
    double gainPerRound = 0.0;

    for (count l = 0; l < maxIter; l++) {
        auto start = std::chrono::high_resolution_clock::now();
        gainPerRound = 0.0;
        graph.parallelForNodes([&](node u) {
            if (!vaff[u])
                return;
            vaff[u] = false;
            auto [bestCommunity, gain] =
                getBestCommunity(graph, u, rho_total, communityMemberships, communitySizes,
                                 edgeCommunityMemberships, edgeCommunityVolumes);

            if (bestCommunity != communityMemberships[u]) {
                updateMemberships(graph, u, bestCommunity, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);

                graph.forNeighborsOf(u, [&](node w) { vaff[w] = true; });

#pragma omp atomic
                gainPerRound += gain;
            }
        });

        auto end = std::chrono::high_resolution_clock::now();
        double duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        INFO("Greedy move phase time for round ", l, ": ", duration, " seconds");
        INFO("Gain in round ", l, ": ", gainPerRound);
        if (gainPerRound < tolerance) {
            break;
        }
    }
};

void HyperLeiden::refineDisconnected(const Hypergraph &graph,
                                     std::vector<count> &communityMemberships,
                                     std::vector<count> &communitySizes,
                                     std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                                     std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes,
                                     std::vector<count> &referenceCommunityMemberships) {

    count maxIter = 100;
    count rho_total = graph.edgeVolume();

    std::vector<bool> vaff(graph.numberOfNodes(), true);
    for (count l = 0; l < maxIter; l++) {
        graph.parallelForNodes([&](node u) {
            if (!vaff[u])
                return;
            if (communitySizes[communityMemberships[u]] != 1) // Only perform for isolated nodes
                return;
            vaff[u] = false;
            auto start = std::chrono::high_resolution_clock::now();
            auto [bestCommunity, gain] = getBestCommunity(
                graph, u, rho_total, communityMemberships, communitySizes, edgeCommunityMemberships,
                edgeCommunityVolumes, referenceCommunityMemberships);
            auto end = std::chrono::high_resolution_clock::now();
            double duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("get best community time for node ", u, ": ", duration, " seconds");
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
                                            std::vector<count> &communityMemberships,
                                            std::vector<count> &communitySizes) {

    // Renumber communities based on prefix sum
    renumberCommunities(communityMemberships, communitySizes);

    // Create new hypergraph with number of nodes = number of communities
    Hypergraph aggHypergraph(communitySizes[communitySizes.size() - 1], 0, true);

    // TODO: maybe there is something more efficient here
    // Iterate over edges of the original hypergraph
    // For each edge, if singleton, add to new hypergraph and track the community id
    // (vector of comm ids) For each edge, if not singleton, create hash of sorted
    // membership vectors (hash combine) and add to a hash database, which maps
    // set-hashes to edge ids.
    // For each hash in the hash database, create a new edge in the new hypergraph and add
    // the corresponding node ids with their combined weights.
    auto convertBitsToNodes = [&](const std::vector<bool> &communityBits) {
        std::vector<node> nodesInEdge;
        for (size_t i = 0; i < communityBits.size(); ++i) {
            if (communityBits[i]) {
                nodesInEdge.push_back(i);
            }
        }
        return nodesInEdge;
    };

    std::unordered_map<std::vector<bool>, std::vector<edgeid>> hashDatabase;
    graph.forEdges([&](edgeid eId) {
        auto nodes = graph.nodesOf(eId);
        bool isSingletonEdge = true;
        if (graph.order(eId) == 0) {
            return; // Skip empty edges
        }
        count indicatorCommunity = nodes.begin()->first; // Get the first node's community id
        nodeweight indicatorWeight = 1;
        std::vector<bool> communityBits(communitySizes[communitySizes.size() - 1], false);

        for (const auto &node : nodes) {
            if (isSingletonEdge && communityMemberships[node.first] != indicatorCommunity) {
                isSingletonEdge = false; // If any node has a different community id, not singleton
            }
            communityBits[communityMemberships[node.first]] = true; // Mark the community bit
            indicatorWeight += node.second;                         // Sum the weights of the nodes
        }

        if (isSingletonEdge) {
            // If singleton edge, add to new hypergraph
            aggHypergraph.addEdge({indicatorCommunity}, false, {indicatorWeight});
        } else {
            // If not singleton, create a hash of the community bits
            auto it = hashDatabase.find(communityBits);
            if (it == hashDatabase.end()) {
                hashDatabase[communityBits] = {eId};
            } else {
                // If found, add edge to the entry
                hashDatabase[communityBits].push_back(eId);
            }
        }

        for (const auto &entry : hashDatabase) {
            // Add edge to hypergraph
            edgeid newEdge = aggHypergraph.addEdge(convertBitsToNodes(entry.first));
            // Gather weights from the edges in the hash database
            for (const auto edgeId : entry.second) {
                auto edgeNodes = graph.nodesOf(edgeId);
                for (const auto &nodeEntry : edgeNodes) {
                    aggHypergraph.updateNodeWeightOf(newEdge, communityMemberships[nodeEntry.first],
                                                     nodeEntry.second);
                }
            }
        }
    });
    return aggHypergraph;
}

// Helper Functions

void HyperLeiden::initializeMemberships(
    const Hypergraph &graph, std::vector<count> &communityMemberships,
    std::vector<count> &communitySizes, std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
    std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) const {
    // In the beginning each node is its own community
    std::iota(std::begin(communityMemberships), std::end(communityMemberships), 0);
    std::fill(std::begin(communitySizes), std::end(communitySizes), 1);

    graph.parallelForEdges([&](edgeid eId, edgeweight /*/*/) {
        auto handleMemberships = edgeCommunityMemberships[eId].makeHandle();
        auto handleVolumes = edgeCommunityVolumes[eId].makeHandle();
        auto nodes = graph.nodesOf(eId);
        for (auto &it : nodes) {
            handleMemberships->insert(communityMemberships[it.first], 1);
            handleVolumes->insert(communityMemberships[it.first], it.second);
        }
    });
}

// TODO: Maybe doing this directly is more efficient than saving the infos in communities variable
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
HyperLeiden::getBestCommunity(const Hypergraph &graph, node v, count rho_total,
                              const std::vector<count> &communityMemberships,
                              const std::vector<count> &communitySizes,
                              std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                              std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes,
                              const std::vector<count> &referenceCommunityMemberships) const {
    // auto [oldCommunity, oldSize] =
    //     std::make_pair(communityMemberships[v], communitySizes[communityMemberships[v]]);
    auto oldCommunity = communityMemberships[v];
    auto oldSize = communitySizes[oldCommunity];
    auto bestCommunity = oldCommunity;
    auto bestGain = std::numeric_limits<double>::lowest();
    // auto [bestCommunity, bestGain] = std::make_pair(0, 0.0);
    // auto neighboringCommunities =
    //     gatherNeighboringCommunities(graph, v, communityMemberships, communitySizes);
    // for (const auto &[c, size] : neighboringCommunities) {
    //     double gain = deltaHCPM(graph, v, oldCommunity, c, oldSize, size, communityMemberships,
    //                             communitySizes, edgeCommunityMemberships, edgeCommunityVolumes);
    //     if (gain > bestGain) {
    //         bestGain = gain;
    //         bestCommunity = c;
    //     }
    // }
    // Measures to be faster doing this directly
    graph.forNeighborsOf(v, [&](node w) {
        if (communityMemberships[v] != communityMemberships[w]) {
            double gain =
                deltaHCPM(graph, v, rho_total, oldCommunity, communityMemberships[w], oldSize,
                          communitySizes[communityMemberships[w]], communityMemberships,
                          communitySizes, edgeCommunityMemberships, edgeCommunityVolumes);
            if (gain > bestGain) {
                bestGain = gain;
                bestCommunity = communityMemberships[w];
            }
        }
    });

    // TODO: maybe default gain not 0.0
    return {bestCommunity, bestGain};
}

double HyperLeiden::deltaHCPM(const Hypergraph &graph, node v, count rho_total, count c1, count c2,
                              count c1Size, count c2Size,
                              const std::vector<count> &communityMemberships,
                              const std::vector<count> &communitySizes,
                              std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                              std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) const {
    double delta_n = 0.0;
    double delta_e = 0.0;
    count iFactor = 0.0;
    count jFactor = 0.0;
    count c1eSize = 0.0;
    count c2eSize = 0.0;

    if (c1 != c2) {
        // delta_n dominates delta_e? Certainly if rho_total >> eWeight
        // If so, we can use this information to speed up the computation
        // delta_n = gamma * rho_total * ((1 << c2Size) - (1 << (c1Size - 1)));
        delta_n = gamma * ((1 << c2Size) - (1 << (c1Size - 1)));
        // For edges of should be faster, current code takes 66 seconds for 1 greedy move iteration
        // graph.forEdges(
        graph.forEdgesOf(v, [&](edgeid eId, edgeweight eWeight) {
            if (graph.hasNode(v, eId)) {
                nodeweight nWeight = graph.getNodeWeightOf(v, eId);

                auto handleMemberships = edgeCommunityMemberships[eId].makeHandle();
                auto handleVolumes = edgeCommunityVolumes[eId].makeHandle();

                jFactor = handleVolumes->find(c2) + 2 * nWeight;
                iFactor = handleVolumes->find(c1) + nWeight;
                handleVolumes->find(c2) == Aux::ParallelHashMap::ht_invalid_value
                    ? jFactor = 2 *nWeight
                    : jFactor = handleVolumes->find(c2) + 2 * nWeight;
                handleVolumes->find(c1) == Aux::ParallelHashMap::ht_invalid_value
                    ? iFactor = nWeight
                    : iFactor = handleVolumes->find(c1) + nWeight;
                handleMemberships->find(c2) == Aux::ParallelHashMap::ht_invalid_value
                    ? c2eSize = 0
                    : c2eSize = handleMemberships->find(c2);
                // c1 should be present in the edge. Otherwise v could not be part of it
                handleMemberships->find(c1) == Aux::ParallelHashMap::ht_invalid_value
                    ? c1eSize = 0
                    : c1eSize = handleMemberships->find(c1);
                // c2eSize = handleMemberships->find(c2);
                // c1eSize = handleMemberships->find(c1);

                // This formula looks fishy ... why c1esize - 1?
                delta_e += eWeight * ((1 << c2eSize) * jFactor - (1 << (c1eSize - 1)) * iFactor);
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

void HyperLeiden::renumberCommunities(std::vector<count> &communityMemberships,
                                      std::vector<count> &communitySizes) {
    // Renumber communities based on prefix sum. This effectively destroys the old
    // information in communitySizes.
    count numCommunities = 0;
    for (size_t i = 0; i < communitySizes.size(); i++) {
        count entry = communitySizes[i];
        communitySizes[i] = numCommunities;
        if (entry > 0) {
            numCommunities++;
        }
    }

    // Update communityMemberships to reflect the new community ids. This effectively
    // destroy the old information in communityMemberships.
    for (size_t i = 0; i < communityMemberships.size(); i++) {
        communityMemberships[i] = communitySizes[communityMemberships[i]];
    }
};

std::vector<node> HyperLeiden::createMapping(const std::vector<count> &communityMemberships,
                                             const std::vector<count> &communitySizes) const {
    std::vector<node> mapping(communitySizes[communitySizes.size() - 1]);
    for (size_t i = 0; i < communityMemberships.size(); i++) {
        mapping[communityMemberships[i]] = i;
    }
    return mapping;
}

void HyperLeiden::flattenPartition() {
    if (mappings.empty()) {
        return;
    }
    // Create a new partition with size |V(G)| (the fine/bigger Graph)
    Partition flattenedPartition(G->numberOfNodes());
    flattenedPartition.setUpperBound(G->numberOfNodes());
    int i = mappings.size() - 1;
    std::vector<node> &lower = mappings[i--];
    while (i >= 0) {
        // iteratively "resolve" (i.e compose) mappings. Let "lower" be a mapping thats
        // below "higher" in the hierarchy (i.e. of a later aggregation)
        // If higher[index] = z and lower[z] = x then set higher[index] = x
        std::vector<node> &upper = mappings[i--];
        for (auto &idx : upper) {
            idx = lower[idx];
        }
        lower = upper;
    }
    G->parallelForNodes([&](node a) { flattenedPartition[a] = lower[a]; });
    flattenedPartition.compact(true);
    result = flattenedPartition;
    mappings.clear();
}
} // namespace NetworKit
