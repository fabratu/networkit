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

void HyperLeiden::run(){
    // nothing here yet
};

void HyperLeiden::greedyMovePhase(const Hypergraph &graph, std::vector<count> &communityMemberships,
                                  Aux::HTCustodian &communitySizes,
                                  Aux::HTCustodian &edgeCommunityMemberships){
    // nothing here yet
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
} // namespace NetworKit
