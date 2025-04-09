/*
 * HypergraphLeiden.hpp
 *
 * Created on: 07.2024
 *     Author: Isaline PLAID
 */

#ifndef NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_
#define NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_

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
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/structures/Partition.hpp>

// Hyperedge weight function depending on their order
double weight_1(double x) {
    return x / 2.0;
}

namespace NetworKit {

class HypergraphLeiden final : public CommunityDetectionAlgorithm<Hypergraph> {
public:
    /**
     *
     * @param graph A networkit hypergraph
     * @param iterations Number of Leiden Iterations to be run
     * @param randomize Randomize node order?
     * @param gamma Resolution parameter
     * @param gamma_cut Resolution parameter for refinement phase, by default = 1 / vol(graph)
     * @param type_contribution Type of modularity chosen : 00=strict , or 01=majority , or
     * 10=strict/partially weighted , or 11=majority/partially weighted.
     * @param weightFun Hyperedge weight function depending on their order
     */
    HypergraphLeiden(const Hypergraph &graph, int iterations = 3, bool randomize = true,
                     double gamma = 1, double gamma_cut = none, int type_contribution = 1,
                     double (*weightFun)(double) = &weight_1);

    void run() override;

    int VECTOR_OVERSIZE = 10000;

private:
    void calculateVolumesHypergraph(const Hypergraph &graph);

    void MoveHypergraph(const Hypergraph &graph, const Partition &zeta);

    double deltaModHypergraph(const Hypergraph &graph, const Partition &zeta, index S,
                              const Partition &p, index c, index target_c);

    Partition RefineHypergraph(const Hypergraph &graph, const Partition &zeta);

    double HypergraphCut(const Hypergraph &graph, const Partition &zeta, index S);

    double GraphVolume; // vol(V)
    double GraphVolume_unweighted;

    std::vector<double> communityVolumes_1; // volume of each community
    std::vector<double> communityVolumes_2;

    std::vector<double> communityVolumes_1_unweighted; // volume of each community
    std::vector<double> communityVolumes_2_unweighted;

    int d; // max size of hyperedge

    std::vector<double> EdgeSizeWeight; // vector of sum of weight of edges of size i (i=0 to d)

    double totalEdgeWeight;

    static constexpr int WORKING_SIZE = 1000;

    double gamma; // Resolution parameter

    double gamma_cut;

    int type_contribution;

    bool changed;

    int numberOfIterations;

    Aux::SignalHandler handler;

    bool random;

    std::vector<std::set<index>> NeighborComm;

    int step = 0.0;

    double (*weightFun)(double);
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_
