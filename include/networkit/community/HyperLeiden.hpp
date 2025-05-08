/*
 * HyperLeiden.hpp
 *
 * Created on: 09.04.2025
 *     Author: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_COMMUNITY_HYPER_LEIDEN_HPP_
#define NETWORKIT_COMMUNITY_HYPER_LEIDEN_HPP_

#include <atomic>
#include <condition_variable>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <thread>

#include <tlx/unused.hpp>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/DynamicParallelHashTable.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class HyperLeiden final : public CommunityDetectionAlgorithm<Hypergraph> {
public:
    enum class CoarseningStrategy { STANDARD, MAJORITY_CONTRIBUTION };
    enum class RefinementStrategy { LEIDEN, DISCONNECTED };

    /**
     *
     * @param hGraph A networkit hypergraph
     * @param iterations Number of HyperLeiden ierations to be run
     * @param gamma Resolution parameter
     * @param refinementStrategy Set strategy for the refinement phase. "DISCONNECTED"
     * guarantees connected communities, while "LEIDEN" uses gamma together with the cut computation
     * to split communities after each greedy move phase.
     * @param coarseningStrategy Set strategy for coarsening
     */
    HyperLeiden(const Hypergraph &hGraph, int numberOfIterations = 3, double gamma = 1,
                RefinementStrategy refinementStrategy = RefinementStrategy::DISCONNECTED,
                CoarseningStrategy coarseningStrategy = CoarseningStrategy::STANDARD);

    void run() override;

private:
    // Main Functions
    // maps to "leiderMoveOmpW" in GVE-Leiden
    void greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                         std::vector<count> &communitySizes,
                         Aux::HTCustodian &edgeCommunityMemberships);

    void refineDisconnected(const Hypergraph &graph, std::vector<count> &communityMemberships,
                            std::vector<count> &tmpCommunitySizes,
                            std::vector<count> &referenceCommunityMemberships);

    Hypergraph aggregateHypergraph(const Hypergraph &graph,
                                   const std::vector<count> &communityMemberships);

    // Helper Functions
    void initializeMemberships(std::vector<count> &communityMemberships,
                               std::vector<count> &communitySizes);
    // Check if needed
    void initializeMembershipsWithPartition(std::vector<count> &communityMemberships,
                                            std::vector<count> &communitySizes,
                                            const Partition &baseClustering);
    // maps to "leidenScanCommunityW" in GVE-Leiden
    std::unordered_map<count, count>
    gatherNeighboringCommunities(const Hypergraph &graph, node v,
                                 const std::vector<count> &communityMemberships,
                                 const std::vector<count> &communitySizes) const;

    // maps to "leidenChooseBestCommunity" in GVE-Leiden
    std::pair<count, double>
    getBestCommunity(const Hypergraph &graph, node v,
                     const std::vector<count> &communityMemberships,
                     const std::vector<count> &communitySizes,
                     const std::vector<count> &referenceCommunityMemberships = {}) const;

    // maps to "deltaModularity" in GVE-Leiden
    double deltaHCPM(count c1, count c2, count c1Size, count c2Size,
                     const std::vector<count> &communityMemberships,
                     const std::vector<count> &communitySizes) const;

    // maps to "leidenChangeCommunity" in GVE-Leiden
    template <bool Refine = false>
    bool updateMemberships(node v, count bestCommunity, std::vector<count> &communityMemberships,
                           std::vector<count> &communitySizes) {
        count currentMembership = communityMemberships[v];
        if (Refine) {
            count tmpCSize = 0;
#pragma omp atomic capture
            {
                tmpCSize = communitySizes[currentMembership];
                communitySizes[currentMembership]--; // TODO: add node weights
            }
            if (tmpCSize > 1) { // TODO: add node weights
#pragma omp atomic
                communitySizes[currentMembership]++;
                return false;
            }
        } else {
#pragma omp atomic
            communitySizes[currentMembership]--; // TODO: add node weights
        }
#pragma omp atomic
        communitySizes[bestCommunity]++; // TODO: add node weights
        communityMemberships[v] = bestCommunity;
        return true;
    }

    // Hyperparameter
    double gamma;              // Resolution parameter
    double tolerance = 0.0001; // Tolerance for the greedy move phase
    int numberOfIterations;
    RefinementStrategy refinementStrategy;
    CoarseningStrategy coarseningStrategy;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_
