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
    // Functions
    void greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                         Aux::HTCustodian &communitySizes,
                         Aux::HTCustodian &edgeCommunityMemberships);

    void refineDisconnected(const Hypergraph &graph, std::vector<count> &communityMemberships);

    Hypergraph aggregateHypergraph(const Hypergraph &graph,
                                   const std::vector<count> &communityMemberships);

    // Hyperparameter
    double gamma; // Resolution parameter
    int numberOfIterations;
    RefinementStrategy refinementStrategy;
    CoarseningStrategy coarseningStrategy;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_
