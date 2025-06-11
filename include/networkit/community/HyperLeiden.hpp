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
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/ParallelHashMap.hpp>
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
                double tolerance = 1e-4,
                RefinementStrategy refinementStrategy = RefinementStrategy::DISCONNECTED,
                CoarseningStrategy coarseningStrategy = CoarseningStrategy::STANDARD);

    void run() override;

private:
    // Main Functions
    // maps to "leiderMoveOmpW" in GVE-Leiden
    void greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                         std::vector<count> &communitySizes,
                         std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                         std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes);

    void refineDisconnected(const Hypergraph &graph, std::vector<count> &communityMemberships,
                            std::vector<count> &communitySizes,
                            std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                            std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes,
                            std::vector<count> &referenceCommunityMemberships);

    Hypergraph aggregateHypergraph(const Hypergraph &graph,
                                   std::vector<count> &communityMemberships,
                                   std::vector<count> &communitySizes);

    std::vector<node> createMapping(const std::vector<count> &communityMemberships,
                                    const std::vector<count> &communitySizes) const;

    // Helper Functions
    void initializeMemberships(const Hypergraph &graph, std::vector<count> &communityMemberships,
                               std::vector<count> &communitySizes,
                               std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                               std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) const;

    // maps to "leidenScanCommunityW" in GVE-Leiden
    std::unordered_map<count, count>
    gatherNeighboringCommunities(const Hypergraph &graph, node v,
                                 const std::vector<count> &communityMemberships,
                                 const std::vector<count> &communitySizes) const;

    // maps to "leidenChooseBestCommunity" in GVE-Leiden
    std::pair<count, double>
    getBestCommunity(const Hypergraph &graph, node v, count rho_total,
                     const std::vector<count> &communityMemberships,
                     const std::vector<count> &communitySizes,
                     std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                     std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes,
                     const std::vector<count> &referenceCommunityMemberships = {}) const;

    // maps to "deltaModularity" in GVE-Leiden
    double deltaHCPM(const Hypergraph &graph, node v, count rho_total, count c1, count c2,
                     count c1Size, count c2Size, const std::vector<count> &communityMemberships,
                     const std::vector<count> &communitySizes,
                     std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                     std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) const;

    // maps to "leidenRenumberCommunitiesW" in GVE-Leiden
    void renumberCommunities(std::vector<count> &communityMemberships,
                             std::vector<count> &communitySizes);

    // maps to "leidenChangeCommunity" in GVE-Leiden
    template <bool Refine = false>
    bool updateMemberships(const Hypergraph &graph, node v, count bestCommunity,
                           std::vector<count> &communityMemberships,
                           std::vector<count> &communitySizes,
                           std::vector<Aux::ParallelHashMap> &edgeCommunityMemberships,
                           std::vector<Aux::ParallelHashMap> &edgeCommunityVolumes) const {
        count currentMembership = communityMemberships[v];
        if (Refine) {
            count tmpCSize = 0;
            bool validUpdate = true;
#pragma omp atomic capture
            {
                tmpCSize = communitySizes[currentMembership];
                communitySizes[currentMembership]--;
            }
            if (tmpCSize > 1) {
#pragma omp atomic
                communitySizes[currentMembership]++;
                validUpdate = false;
            }
            if (!validUpdate) {
                return false;
            }
            auto handleMemberships = edgeCommunityMemberships[v].makeHandle();
            auto handleVolumes = edgeCommunityVolumes[v].makeHandle();
            auto currentMemberships = handleMemberships->find(v);
            auto currentVolumes = handleVolumes->find(v);
            handleMemberships->update(currentMembership, currentMemberships - 1);
            handleVolumes->update(currentMembership, currentVolumes - graph.getNodeWeight(v));
        } else {
            auto edgesOf = graph.edgesOf(v);
            for (auto &it : edgesOf) {
                auto handleMemberships = edgeCommunityMemberships[it].makeHandle();
                auto handleVolumes = edgeCommunityVolumes[it].makeHandle();
                auto currentMemberships = handleMemberships->find(v);
                auto currentVolumes = handleVolumes->find(v);
                handleMemberships->update(currentMembership, currentMemberships - 1);
                handleVolumes->update(currentMembership, currentVolumes - graph.getNodeWeight(v));
            }
#pragma omp atomic
            communitySizes[currentMembership]--;
        }
        auto edgesOf = graph.edgesOf(v);
        for (auto &it : edgesOf) {
            auto handleMemberships = edgeCommunityMemberships[it].makeHandle();
            auto handleVolumes = edgeCommunityVolumes[it].makeHandle();
            auto bestMemberships = handleMemberships->find(v);
            auto bestVolumes = handleVolumes->find(v);
            if (v == Aux::ParallelHashMap::ht_invalid_key) {
                handleMemberships->insert(bestCommunity, 1);
                handleVolumes->insert(bestCommunity, graph.getNodeWeight(v));
            } else {
                handleMemberships->update(bestMemberships, bestMemberships + 1);
                handleVolumes->update(bestMemberships, bestVolumes + graph.getNodeWeight(v));
            }
        }
#pragma omp atomic
        communitySizes[bestCommunity]++;
        communityMemberships[v] = bestCommunity;
        return true;
    }

    std::vector<bool> communityExists(const Hypergraph &graph, node v, count bestCommunity,
                                      const std::vector<count> &communityMemberships,
                                      const std::vector<count> &communitySizes) const;

    void flattenPartition();

    // Hyperparameter
    double gamma;     // Resolution parameter
    double tolerance; // Tolerance for the greedy move phase
    int numberOfIterations;
    RefinementStrategy refinementStrategy;
    CoarseningStrategy coarseningStrategy;
    std::vector<std::vector<node>> mappings;
}; // namespace NetworKit

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_
