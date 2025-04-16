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
                                  std::vector<count> &communitySizes,
                                  Aux::HTCustodian &edgeCommunityMemberships) {

    auto deltaHCPM = [&](count c1, count c2, count c1Size, count c2Size) {
        double delta = 0.0;

        if (c1 != c2) {
            delta -= 1.0;
            // delta = (communitySizes[c1] - communitySizes[c2]) / graph.numberOfEdges();
            // delta += (edgeCommunityMemberships[c1] - edgeCommunityMemberships[c2])
            //          / graph.numberOfEdges();
        }
        return delta;
    };

    auto gatherNeighboringCommunities = [&](node v) {
        std::unordered_set<std::pair<count, count>> communities;
        graph.forNeighborsOf(v, [&](node w) {
            if (communityMemberships[v] != communityMemberships[w]) {
                communities.insert(
                    {communityMemberships[w], communitySizes[communityMemberships[w]]});
            }
        });
        return communities;
    };

    auto getBestCommunity = [&](node v) -> std::pair<count, double> {
        count oldCommunity = communityMemberships[v];
        count bestCommunity = communityMemberships[v];
        double maxGain = 0.0;
        auto neighboringCommunities = gatherNeighboringCommunities(v);

        for (auto &community : neighboringCommunities) {
            double gain =
                deltaHCPM(communityMemberships[v], oldCommunity, community.first, community.second);
            if (gain > maxGain) {
                maxGain = gain;
                bestCommunity = community.first;
            }
        }
        return {bestCommunity, maxGain};
    };

    std::vector<bool> vaff(graph.numberOfNodes(), false);
    count numAffected = graph.numberOfNodes();

    do {
        graph.parallelForNodes([&](node u) {
            if (!vaff[u])
                return;
            vaff[u] = false;
            numAffected--;
            auto [bestCommunity, gain] = getBestCommunity(u);
            if (bestCommunity != communityMemberships[u]) {

                // communityMemberships[u] = bestCommunity;
                // communitySizes[bestCommunity]++;
                // communitySizes[communityMemberships[u]]--;
                // edgeCommunityMemberships[bestCommunity]++;
                // edgeCommunityMemberships[communityMemberships[u]]--;
            }
        });
        // #pragma omp parallel for schedule(dynamic, 2048) reduction(+ : el)
        //         for (K u = 0; u < S; ++u) {
        //             int t = omp_get_thread_num();
        //             if (!x.hasVertex(u))
        //                 continue;
        //             if (!fa(u) || !vaff[u])
        //                 continue;
        //             if (REFINE && ctot[vcom[u]] > vtot[u])
        //                 continue;
        //             leidenClearScanW(*vcs[t], *vcout[t]);
        //             leidenScanCommunitiesW<false, REFINE>(*vcs[t], *vcout[t], x, u, vcom, vcob);
        //             auto [c, e] = leidenChooseCommunity(x, u, vcom, vtot, ctot, *vcs[t],
        //             *vcout[t], M, R); if (c && leidenChangeCommunityOmpW<REFINE>(vcom, ctot, x,
        //             u, c, vtot))
        //                 x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); });
        //             vaff[u] = B();
        //             el += e; // l1-norm
        //         }
        //         if (REFINE || fc(el, l++))
        //             break;

    } while (numAffected);
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
