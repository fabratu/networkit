/*
 * HyperLeiden.cpp
 *
 * Created on: 09.04.2025
 *     Author: Fabian Brandt-Tumescheit
 */

#include <boost/functional/hash.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/HyperLeiden.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/structures/Partition.hpp>

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
    auto lastNumberOfNodes = G->numberOfNodes();

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

            // result.reset(G->upperNodeIdBound(), 0);
            // result.setUpperBound(G->upperNodeIdBound());
            // // Partition zeta(communityMemberships.size());
            // // zeta.allToSingletons();
            // G->parallelForNodes([&](node u) { result[u] = communityMemberships[u]; });
            // Modularity mod;
            // double quality = mod.getQualityHypergraph(result, *G);
            // INFO("Modularity quality of initial partition: ", quality);

            // Refine disconnected communities
            start = std::chrono::high_resolution_clock::now();
            if (refinementStrategy == RefinementStrategy::DISCONNECTED) {
                std::vector<count> tmpCommunityMemberships(G->upperNodeIdBound(), 0);
                std::vector<count> tmpCommunitySizes(G->upperNodeIdBound(), 0);
                std::vector<Aux::ParallelHashMap> tmpEdgeCommunityMemberships(
                    G->upperEdgeIdBound());
                std::vector<Aux::ParallelHashMap> tmpEdgeCommunityVolumes(G->upperEdgeIdBound());
                refineDisconnected(*G, tmpCommunityMemberships, tmpCommunitySizes,
                                   tmpEdgeCommunityMemberships, tmpEdgeCommunityVolumes,
                                   communityMemberships);
            }

            end = std::chrono::high_resolution_clock::now();
            duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("Refinement time: ", duration, " seconds");
            // result.reset(G->upperNodeIdBound(), 0);
            // result.setUpperBound(G->upperNodeIdBound());
            // // Partition zeta(communityMemberships.size());
            // // zeta.allToSingletons();
            // for (node u = 0; u < communityMemberships.size(); ++u) {
            //     result.moveToSubset(communityMemberships[u], u);
            // }
            // quality = mod.getQualityHypergraph(result, *G);
            // INFO("Modularity quality of initial partition: ", quality, " with gamma = ", gamma);

            // exit(0);
            // Aggregate hypergraph
            start = std::chrono::high_resolution_clock::now();
            currentG = aggregateHypergraph(*G, communityMemberships, communitySizes);
            end = std::chrono::high_resolution_clock::now();
            duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            INFO("Aggregation time: ", duration, " seconds");
            INFO("Graph after first pass has ", currentG.numberOfNodes(), " nodes and ",
                 currentG.numberOfEdges(), " edges.");

            // Print all hyperedges of the graph
            // INFO("Printing all hyperedges of the graph:");
            // currentG.forEdges([&](edgeid eId) {
            //     auto nodes = currentG.nodesOf(eId);
            //     std::string edgeStr = "Edge " + std::to_string(eId) + ": {";
            //     bool first = true;
            //     for (const auto &nodePair : nodes) {
            //         if (!first)
            //             edgeStr += ", ";
            //         edgeStr += std::to_string(nodePair.first)
            //                    + "(w=" + std::to_string(nodePair.second) + ")";
            //         first = false;
            //     }
            //     edgeStr += "}";
            //     INFO(edgeStr);
            // });
            // Create mapping
            mappings.push_back(createMapping(communityMemberships, communitySizes));
            isFirst = false;
            if (currentG.numberOfNodes() == lastNumberOfNodes || currentG.numberOfNodes() <= 2) {
                INFO("Stopping early, no further aggregation possible.");
                break;
            }
            lastNumberOfNodes = currentG.numberOfNodes();
        } else {
            INFO("Starting iteration ", pass, " of HyperLeiden");
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
            INFO("Graph after current pass has ", currentG.numberOfNodes(), " nodes and ",
                 currentG.numberOfEdges(), " edges.");

            // Print all hyperedges of the graph
            // INFO("Printing all hyperedges of the graph:");
            // currentG.forEdges([&](edgeid eId) {
            //     auto nodes = currentG.nodesOf(eId);
            //     std::string edgeStr = "Edge " + std::to_string(eId) + ": {";
            //     bool first = true;
            //     for (const auto &nodePair : nodes) {
            //         if (!first)
            //             edgeStr += ", ";
            //         edgeStr += std::to_string(nodePair.first)
            //                    + "(w=" + std::to_string(nodePair.second) + ")";
            //         first = false;
            //     }
            //     edgeStr += "}";
            //     INFO(edgeStr);
            // });

            // Create mapping
            mappings.push_back(createMapping(communityMemberships, communitySizes));
            INFO("Created mapping for pass ", pass);
            if (currentG.numberOfNodes() == lastNumberOfNodes || currentG.numberOfNodes() <= 4) {
                INFO("Stopping early, no further aggregation possible.");
                break;
            }
            lastNumberOfNodes = currentG.numberOfNodes();
        }
    }

    // TODO: unroll mappings to get final community memberships
    // INFO("Final mappings: ", Aux::toString(mappings));
    flattenPartition();

    hasRun = true;
}

void HyperLeiden::greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                                  std::vector<count> &communitySizes,
                                  std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                                  std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) {

    count maxIter = 1000;
    count rho_total = graph.edgeVolume();

    std::vector<bool> vaff(graph.numberOfNodes(), true);
    double gainPerRound = 0;
    double lastGain = 0.0;

    for (count l = 0; l < maxIter; l++) {
        auto start = std::chrono::high_resolution_clock::now();
        gainPerRound = 0.0;

        graph.parallelForNodesInRandomOrder([&](node u) {
            if (!vaff[u])
                return;
            vaff[u] = false;
            auto [bestCommunity, gain] =
                getBestCommunity(graph, u, rho_total, communityMemberships, communitySizes,
                                 edgeCommunityMemberships, edgeCommunityVolumes);

            if (bestCommunity != communityMemberships[u] && gain > 0) {
                // INFO("Node ", u, " moves to community ", bestCommunity, " with gain ", gain,
                //      " in round ", l);
                updateMemberships(graph, u, bestCommunity, communityMemberships, communitySizes,
                                  edgeCommunityMemberships, edgeCommunityVolumes);

                graph.forNeighborsOf(u, [&](node w) { vaff[w] = true; });

                // if (gain > 100000) {
                //     int64_t tempGain = gain;
                //     count tempBestCommunity = bestCommunity;
                //     INFO("Node ", u, " moved to community ", bestCommunity, " with gain ", gain,
                //          " in round ", l);
                //     exit(1);
                //     l = 100;
                // }
#pragma omp atomic
                gainPerRound += gain;
                // INFO("Node ", u, " moved to community ", bestCommunity, " with gain ", gain);
            }
        });

        auto end = std::chrono::high_resolution_clock::now();
        double duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        INFO("Greedy move phase time for round ", l, ": ", duration, " seconds");
        INFO("Gain in round ", l, ": ", gainPerRound);
        if (gainPerRound < tolerance || gainPerRound <= lastGain) {
            INFO("Stopping greedy move phase");
            break;
        }
        lastGain = gainPerRound;
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
        graph.parallelForNodesInRandomOrder([&](node u) {
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
            // INFO("get best community time for node ", u, ": ", duration, " seconds");
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
    count numCommunities = renumberCommunities(communityMemberships, communitySizes);
    // INFO("Renumbered communities: ", Aux::toString(communityMemberships));
    // INFO("Community sizes: ", Aux::toString(communitySizes));
    INFO("Number of communities after renumbering: ", numCommunities);

    // Create new hypergraph with number of nodes = number of communities
    Hypergraph aggHypergraph(numCommunities, 0, true);

    struct CustomHasher {
        // noexcept is recommended, but not required
        std::size_t operator()(const std::vector<bool> &s) const /*noexcept*/
        {
            return boost::hash_range(s.begin(), s.end());
        }
    };

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

    // std::unordered_map<std::vector<bool>, std::vector<edgeid>> hashDatabase;
    std::unordered_map<std::vector<bool>, std::vector<edgeid>, CustomHasher> hashDatabase;
    std::vector<edgeid> singletonDatabase(numCommunities, std::numeric_limits<edgeid>::max());

    INFO("Starting aggregation of hypergraph with ", aggHypergraph.numberOfNodes(), " nodes and ",
         aggHypergraph.numberOfEdges(), " edges.");

    graph.forEdges([&](edgeid eId) {
        if (graph.order(eId) == 0) {
            return; // Skip empty edges
        }

        auto nodes = graph.nodesOf(eId);
        bool isSingletonEdge = true;

        count indicatorCommunity =
            communityMemberships[nodes.begin()->first]; // Get the first node's community id
        nodeweight indicatorWeight = 0.0;
        std::vector<bool> communityBits(numCommunities, false);

        for (const auto &node : nodes) {
            if (isSingletonEdge && communityMemberships[node.first] != indicatorCommunity) {
                isSingletonEdge = false; // If any node has a different community id, not singleton
            }
            communityBits[communityMemberships[node.first]] = true; // Mark the community bit
            indicatorWeight += node.second;                         // Sum the weights of the nodes
        }

        if (isSingletonEdge) {
            // If singleton edge, add to new hypergraph or update weight
            if (singletonDatabase[indicatorCommunity] == std::numeric_limits<edgeid>::max()) {

                // If this is the first singleton edge for this community, add it
                singletonDatabase[indicatorCommunity] =
                    aggHypergraph.addEdge({indicatorCommunity}, false, {indicatorWeight});
                // INFO("Adding singleton edge with community ", indicatorCommunity, " and weight ",
                //      indicatorWeight, " for edge id: ", eId,
                //      " new edge id: ", singletonDatabase[indicatorCommunity]);

            } else {
                // INFO("Updating singleton edge with community ", indicatorCommunity, " and weight
                // ",
                //      indicatorWeight, " for edge id: ", eId,
                //      " new edge id: ", singletonDatabase[indicatorCommunity]);
                // If already exists, update the weight
                aggHypergraph.updateNodeWeightOf(
                    indicatorCommunity, singletonDatabase[indicatorCommunity], indicatorWeight);
            }
        } else {
            // If not singleton, create a hash of the community bits
            auto it = hashDatabase.find(communityBits);
            if (it == hashDatabase.end()) {
                // auto communityHash = std::hash<std::vector<bool>>{}(communityBits);
                // auto communityHash2 = boost::hash_range(
                //     communityBits.begin(), communityBits.end()); // Combine bits into a hash
                // INFO("Creating new entry in hash database for community bits: ",
                //      Aux::toString(communityBits), " with hash: ", communityHash, " and ",
                //      communityHash2, " for edge id: ", eId);
                hashDatabase[communityBits] = {eId};
            } else {
                // If found, add edge to the entry
                hashDatabase[communityBits].push_back(eId);
            }
        }
    });
    INFO("After singleton the graph has ", aggHypergraph.numberOfNodes(), " nodes and ",
         aggHypergraph.numberOfEdges(), " edges.");
    // INFO("Singleton edges in the database: ", Aux::toString(singletonDatabase));
    // INFO("Each singleton edge has node id and weight: ");
    // aggHypergraph.forEdges([&](edgeid eId) {
    //     auto nodes = aggHypergraph.nodesOf(eId);
    //     std::string edgeStr = "Edge " + std::to_string(eId) + ": {";
    //     bool first = true;
    //     for (const auto &nodeEntry : nodes) {
    //         if (!first)
    //             edgeStr += ", ";
    //         edgeStr +=
    //             std::to_string(nodeEntry.first) + "(w=" + std::to_string(nodeEntry.second) + ")";
    //         first = false;
    //     }
    //     edgeStr += "}";
    //     INFO(edgeStr);
    // });
    INFO("Hash database has ", hashDatabase.size(), " entries.");
    for (const auto &entry : hashDatabase) {
        // INFO("Processing hash entry with community bits: ", Aux::toString(entry.first),
        //      " and edges: ", Aux::toString(entry.second));
        // Add edge to hypergraph
        edgeid newEdge = aggHypergraph.addEdge(convertBitsToNodes(entry.first));
        auto edgeNodes = aggHypergraph.nodesOf(newEdge);
        for (const auto &nodeEntry : edgeNodes) {
            aggHypergraph.setNodeWeightOf(nodeEntry.first, newEdge, 0.0);
        }
        // Gather weights from the edges in the hash database
        for (const auto edgeId : entry.second) {
            auto edgeNodes = graph.nodesOf(edgeId);
            for (const auto &nodeEntry : edgeNodes) {
                aggHypergraph.updateNodeWeightOf(communityMemberships[nodeEntry.first], newEdge,
                                                 nodeEntry.second);
            }
        }
    }
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
            handleMemberships.insert(communityMemberships[it.first], 1);
            handleVolumes.insert(communityMemberships[it.first], it.second);
        }
    });
}

// TODO: Maybe doing this directly is more efficient than saving the infos in communities
// variable
std::unordered_map<count, count>
HyperLeiden::gatherNeighboringCommunities(const Hypergraph &graph, node v,
                                          const std::vector<count> &communityMemberships,
                                          const std::vector<count> &communitySizes) const {
    std::unordered_map<count, count> communities;
    graph.forNeighborsOf(v, [&](node w) {
        if (communityMemberships[v] != communityMemberships[w]) {
            communities[communityMemberships[w]] = communitySizes[communityMemberships[w]];
            // communities.insert({communityMemberships[w],
            // communitySizes[communityMemberships[w]]});
        }
    });
    return communities;
};

// TODO: We can make int64_t smaller
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
    double bestGain = std::numeric_limits<double>::lowest();
    // auto [bestCommunity, bestGain] = std::make_pair(0, 0.0);
    auto neighboringCommunities =
        gatherNeighboringCommunities(graph, v, communityMemberships, communitySizes);
    // INFO("Node ", v, " has ", neighboringCommunities.size(), " neighboring communities");
    for (const auto &[c, size] : neighboringCommunities) {
        // INFO("Node ", v, " has neighboring community ", c, " with size ", size);
    }
    for (const auto &[c, size] : neighboringCommunities) {
        // if (v == 3319 && c == 5739) {
        //     INFO("Envoking deltaHCPM");
        // }
        double gain =
            deltaHCPM(graph, v, rho_total, oldCommunity, c, oldSize, size, communityMemberships,
                      communitySizes, edgeCommunityMemberships, edgeCommunityVolumes);
        // INFO("Node ", v, " gain for community ", c, ": ", gain);
        if (gain > bestGain) {
            bestGain = gain;
            bestCommunity = c;
        }
    }

    // Measures to be faster doing this directly
    // graph.forNeighborsOf(v, [&](node w) {
    //     if (communityMemberships[v] != communityMemberships[w]) {
    //         if (v == 3319 && communityMemberships[w] == 8128) {
    //             INFO("Envoking deltaHCPM");
    //         }
    //         int64_t gain =
    //             deltaHCPM(graph, v, rho_total, oldCommunity, communityMemberships[w],
    //             oldSize,
    //                       communitySizes[communityMemberships[w]], communityMemberships,
    //                       communitySizes, edgeCommunityMemberships, edgeCommunityVolumes);
    //         if (gain > bestGain) {
    //             bestGain = gain;
    //             bestCommunity = communityMemberships[w];
    //         }
    //     }
    // });

    // TODO: maybe default gain not 0.0
    return {bestCommunity, bestGain};
}

double HyperLeiden::deltaHCPM(const Hypergraph &graph, node v, count rho_total, count c1, count c2,
                              count c1Size, count c2Size,
                              const std::vector<count> &communityMemberships,
                              const std::vector<count> &communitySizes,
                              std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                              std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) const {
    int64_t delta_n = 0;
    double delta_e = 0;
    count numEdges = graph.numberOfEdges();

    if (c1 != c2) {
        delta_n = gamma * ((1 << c2Size) - (1 << (c1Size - 1)));
        // INFO("Node ", v, " delta_n: ", delta_n, " for communities ", c1, " and ", c2);

        count c1eSize = 0;
        count c2eSize = 0;
        count membershipValueC1 = 0;
        count membershipValueC2 = 0;
        count volumeValueC1 = 0;
        count volumeValueC2 = 0;
        int64_t left = 0;
        int64_t right = 0;
        nodeweight nWeight = 0.0;

        graph.forEdgesOf(v, [&](edgeid eId, edgeweight eWeight) {
            auto handleMemberships = edgeCommunityMemberships[eId].makeHandle();
            auto handleVolumes = edgeCommunityVolumes[eId].makeHandle();

            membershipValueC1 = handleMemberships.find(c1);
            membershipValueC2 = handleMemberships.find(c2);
            volumeValueC1 = handleVolumes.find(c1);
            volumeValueC2 = handleVolumes.find(c2);
            nWeight = graph.getNodeWeightOf(v, eId);
            left = 0;
            right = 0;

            membershipValueC1 == Aux::ParallelHashMap::ht_invalid_value
                ? c1eSize = 0
                : c1eSize = membershipValueC1;
            membershipValueC2 == Aux::ParallelHashMap::ht_invalid_value
                ? c2eSize = 0
                : c2eSize = membershipValueC2;

            if (c1eSize == 1) {
                right = volumeValueC1;
            } else if (c1eSize > 1) {
                right = (1 << (c1eSize - 2)) * (volumeValueC1 + nWeight);
            }

            if (c2eSize == 0) {
                left = nWeight;
            } else if (c2eSize > 0) {
                left = (1 << (c2eSize - 1)) * (volumeValueC2 + 2 * nWeight);
            }
            // volumeValueC1 != Aux::ParallelHashMap::ht_invalid_value || c1eSize == 1
            //     ? right = volumeValueC1
            //     : right = (1 << (c1eSize - 2)) * (volumeValueC1 + nWeight);
            // volumeValueC2 == Aux::ParallelHashMap::ht_invalid_value
            //     ? left = nWeight
            //     : left = (1 << (c2eSize - 1)) * (volumeValueC2 + 2 * nWeight);

            // int64_t left = static_cast<int64_t>(1 << c2eSize);
            // int64_t right = static_cast<int64_t>(1 << (c1eSize - 1));

            double eVolume = static_cast<double>(graph.edgeVolume(eId));
            // double eFactor = eVolume * (left - right);
            double eFactor = (left - right);

            delta_e += eFactor;
            // INFO("Node ", v, " edge ", eId, " eFactor: ", eFactor, " for communities ", c1, "
            // and
            // ",
            //      c2, " left: ", left, " right: ", right, " with edge weight ", eWeight,
            //      " and eVolume ", eVolume, " c1eSize: ", c1eSize, " c2eSize: ", c2eSize,
            //      " membershipValueC1: ", membershipValueC1,
            //      " membershipValueC2: ", membershipValueC2, " volumeValueC1: ",
            //      volumeValueC1, " volumeValueC2: ", volumeValueC2, " nWeight: ", nWeight, "
            //      delta_e: ", delta_e);
        });

        // c1 should be present in the edge. Otherwise v could not be part of it
        // graph.forEdgesOf(v, [&](edgeid eId, edgeweight eWeight) {
        //     nodeweight nWeight = graph.getNodeWeightOf(v, eId);

        //     auto handleMemberships = edgeCommunityMemberships[eId].makeHandle();
        //     auto handleVolumes = edgeCommunityVolumes[eId].makeHandle();

        //     membershipValueC1 = handleMemberships->find(c1);
        //     membershipValueC2 = handleMemberships->find(c2);
        //     volumeValueC1 = handleVolumes->find(c1);
        //     volumeValueC2 = handleVolumes->find(c2);

        //     membershipValueC1 == Aux::ParallelHashMap::ht_invalid_value
        //         ? c1eSize = 0
        //         : c1eSize = membershipValueC1;
        //     membershipValueC2 == Aux::ParallelHashMap::ht_invalid_value
        //         ? c2eSize = 0
        //         : c2eSize = membershipValueC2;

        //     volumeValueC1 == Aux::ParallelHashMap::ht_invalid_value
        //         ? iFactor = nWeight
        //         : iFactor = volumeValueC1 + nWeight;
        //     volumeValueC2 == Aux::ParallelHashMap::ht_invalid_value
        //         ? jFactor = 2 *nWeight
        //         : jFactor = volumeValueC2 + 2 * nWeight;

        //     int64_t left = static_cast<int64_t>((1 << c2eSize) * jFactor);
        //     int64_t right = static_cast<int64_t>((1 << (c1eSize - 1)) * iFactor);

        //     delta_e += eWeight * (left - right);
        // });
    }
    // INFO("Node ", v, " delta_e: ", delta_e, " for communities ", c1, " and ", c2);
    return delta_e / static_cast<double>(numEdges) - delta_n;
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

count HyperLeiden::renumberCommunities(std::vector<count> &communityMemberships,
                                       std::vector<count> &communitySizes) {
    // Renumber communities based on prefix sum. This effectively destroys the old
    // information in communitySizes.
    // INFO("Community memberships before renumbering: ", Aux::toString(communityMemberships));
    // INFO("Community sizes before renumbering: ", Aux::toString(communitySizes));
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
    return numCommunities;
};

std::vector<node> HyperLeiden::createMapping(const std::vector<count> &communityMemberships,
                                             const std::vector<count> &communitySizes) const {
    std::vector<node> mapping(communityMemberships.size(), 0);
    for (size_t i = 0; i < communityMemberships.size(); i++) {
        mapping[i] = communityMemberships[i];
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
