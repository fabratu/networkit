/*
 * HypergraphCommunityGTest.cpp
 *
 *  Created on: 07.2024
 *      Author: Isaline Plaid
 *              Fabian Brandt-Tumescheit
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/community/Modularity.hpp>
#include <networkit/community/HypergraphLeiden.hpp>
#include <networkit/community/HypergraphLouvain.hpp>

namespace NetworKit {

class HypergraphCommunityGTest : public testing::Test {};

TEST_F(HypergraphCommunityGTest, testLeiden) {
    // Hypergraph creation : 
    bool weighted = true; 
    Hypergraph hg(6, 0, weighted);
    std::vector<node> edge1 = {0, 1, 2};
    edgeid newEdge = hg.addEdge(edge1, true);
    hg.setEdgeWeight(newEdge, 1.5);
    std::vector<node> edge2 = {0, 2};
    newEdge = hg.addEdge(edge2, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge3 = {1, 2};
    newEdge = hg.addEdge(edge3, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge4 = {2, 3};
    newEdge = hg.addEdge(edge4, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge5 = {3, 5};
    newEdge = hg.addEdge(edge5, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge6 = {4, 5};
    newEdge = hg.addEdge(edge6, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge7 = {3, 4, 5};
    newEdge = hg.addEdge(edge7, true);
    hg.setEdgeWeight(newEdge, 1.5);


    Partition q(hg.numberOfNodes());
    q.allToSingletons();
    q.mergeSubsets(q[0], q[2]);
    q.mergeSubsets(q[4], q[5]);
    Partition t(hg.numberOfNodes());
    t.allToSingletons();
    t.mergeSubsets(t[0], t[2]);
    t.mergeSubsets(t[0], t[1]);
    t.mergeSubsets(t[4], t[5]);
    t.mergeSubsets(t[3], t[5]);
    Modularity modularityHypergraph;
    double mod_q = modularityHypergraph.getQualityHypergraph(q, hg, 1);
    double mod_t = modularityHypergraph.getQualityHypergraph(t, hg, 1);
    EXPECT_LE(mod_t, mod_q);

    HypergraphLeiden pl(hg, 3, false);
    pl.run();
    Partition zeta = pl.getPartition();
    double mod_zeta = modularityHypergraph.getQualityHypergraph(zeta, hg, 1);

    ASSERT_TRUE(zeta[0]==2);
    ASSERT_TRUE(zeta[1]==2);
    ASSERT_TRUE(zeta[2]==2);
    ASSERT_TRUE(zeta[3]==5);
    ASSERT_TRUE(zeta[4]==5);
    ASSERT_TRUE(zeta[5]==5);
}

TEST_F(HypergraphCommunityGTest, testLouvain) {
    // Hypergraph creation : 
    bool weighted = true; 
    Hypergraph hg(6, 0, weighted);
    std::vector<node> edge1 = {0, 1, 2};
    edgeid newEdge = hg.addEdge(edge1, true);
    hg.setEdgeWeight(newEdge, 1.5);
    std::vector<node> edge2 = {0, 2};
    newEdge = hg.addEdge(edge2, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge3 = {1, 2};
    newEdge = hg.addEdge(edge3, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge4 = {2, 3};
    newEdge = hg.addEdge(edge4, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge5 = {3, 5};
    newEdge = hg.addEdge(edge5, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge6 = {4, 5};
    newEdge = hg.addEdge(edge6, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge7 = {3, 4, 5};
    newEdge = hg.addEdge(edge7, true);
    hg.setEdgeWeight(newEdge, 1.5);


    Partition q(hg.numberOfNodes());
    q.allToSingletons();
    q.mergeSubsets(q[0], q[2]);
    q.mergeSubsets(q[4], q[5]);
    Partition t(hg.numberOfNodes());
    t.allToSingletons();
    t.mergeSubsets(t[0], t[2]);
    t.mergeSubsets(t[0], t[1]);
    t.mergeSubsets(t[4], t[5]);
    t.mergeSubsets(t[3], t[5]);
    Modularity modularityHypergraph;
    double mod_q = modularityHypergraph.getQualityHypergraph(q, hg, 1);
    double mod_t = modularityHypergraph.getQualityHypergraph(t, hg, 1);
    EXPECT_LE(mod_t, mod_q);

    HypergraphLouvain plouvain(hg, 3, false);
    plouvain.run();
    Partition zeta_bis = plouvain.getPartition();

    ASSERT_TRUE(zeta_bis[0]==2);
    ASSERT_TRUE(zeta_bis[1]==2);
    ASSERT_TRUE(zeta_bis[2]==2);
    ASSERT_TRUE(zeta_bis[3]==5);
    ASSERT_TRUE(zeta_bis[4]==5);
    ASSERT_TRUE(zeta_bis[5]==5);
}

TEST_F(HypergraphCommunityGTest, testModularity) {
    // Hypergraph creation : 
    bool weighted = true; 
    Hypergraph hg(5, 3, weighted);
    std::vector<node> edge1 = {0, 1, 2};
    edgeid newEdge = hg.addEdge(edge1, true);
    hg.setEdgeWeight(newEdge, 4.0);
    std::vector<node> edge2 = {0, 1, 2, 3};
    newEdge = hg.addEdge(edge2, true);
    hg.setEdgeWeight(newEdge, 1.0);
    std::vector<node> edge3 = {1, 4};
    newEdge = hg.addEdge(edge3, true);
    hg.setEdgeWeight(newEdge, 5.0);

    // We partition our hypergraph
    Partition p(hg.numberOfNodes());
    p.allToSingletons();
    p.mergeSubsets(p[0], p[1]);
    p.mergeSubsets(p[1], p[2]);
    p.mergeSubsets(p[1], p[4]);

    // Test modularity value for strict edge contribution
    Modularity modularityHypergraph;
    double mod_0 = modularityHypergraph.getQualityHypergraph(p, hg,1.0,0);
    DEBUG("modularity: ", mod_0);
    ASSERT_TRUE(mod_0 == 0.9 - 2065805. / 2284880.0);

    // Test modularity value for majority edge contribution
    double mod_1 = modularityHypergraph.getQualityHypergraph(p, hg, 1.0,1);
    DEBUG("modularity: ", mod_1);
    ASSERT_TRUE(mod_1 == 1- 2198505.0/2284880.0);

    // Test modularity value for strict edge contribution / partially weighted
    double mod_10 = modularityHypergraph.getQualityHypergraph(p, hg,1.0,10);
    DEBUG("modularity: ", mod_10);
    EXPECT_LE(mod_10 - 0.0000000000000001, 0.9 - 4889.0 / 6561.0); //C++ double should have a floating-point precision of up to 15 digits
    EXPECT_GE(mod_10 + 0.0000000000000001, 0.9 - 4889.0 / 6561.0);

    // Test modularity value for majority edge contribution / partially weighted
    double mod_11 = modularityHypergraph.getQualityHypergraph(p, hg, 1.0,11);
    DEBUG("modularity: ", mod_11);
    ASSERT_TRUE(mod_11 == 1- 9791.0 /10935.0);


    //Testing the value of the modularity gain from moving nodes 0 and 2 into the community of 3
    // For strict edge contribution
    std::set<node> S({0,2});
    double gain_0 = modularityHypergraph.deltaModularityHypergraph(p, hg, S, p[3], 1.0, 0);
    DEBUG("modularity gain: ", gain_0);
    ASSERT_TRUE(gain_0 == -0.4 + 30093.0 / 57122.0);
    // For majority edge contribution
    double gain_1 = modularityHypergraph.deltaModularityHypergraph(p, hg, S, p[3], 1.0, 1);
    DEBUG("modularity gain: ", gain_1);
    ASSERT_TRUE(gain_1 ==0 +13825.0 /57122.0);

    // For strict edge contribution / partially weighted
    double gain_10 = modularityHypergraph.deltaModularityHypergraph(p, hg, S, p[3], 1.0, 10);
    DEBUG("modularity gain: ", gain_10);
    // For majority edge contribution / partially weighted
    double gain_11 = modularityHypergraph.deltaModularityHypergraph(p, hg, S, p[3], 1.0, 11);
    DEBUG("modularity gain: ", gain_11);


    //check the value of gain
    Partition q(hg.numberOfNodes());
    q.allToSingletons();
    q.mergeSubsets(q[0], q[3]);
    q.mergeSubsets(q[0], q[2]);
    q.mergeSubsets(q[1], q[4]);
    double mod_prime_1 = modularityHypergraph.getQualityHypergraph(q, hg, 1.0, 1);
    ASSERT_TRUE(mod_prime_1 == mod_1 + gain_1);
    double mod_prime_0 = modularityHypergraph.getQualityHypergraph(q, hg, 1.0, 0);
    ASSERT_TRUE(mod_prime_0 == mod_0 + gain_0);
    double mod_prime_11 = modularityHypergraph.getQualityHypergraph(q, hg, 1.0, 11);
    EXPECT_LE(mod_prime_11 - 0.0000000000000001, mod_11 + gain_11); //C++ double should have a floating-point precision of up to 15 digits
    EXPECT_GE(mod_prime_11 + 0.0000000000000001, mod_11 + gain_11);
    double mod_prime_10 = modularityHypergraph.getQualityHypergraph(q, hg, 1.0, 10);
    EXPECT_LE(mod_prime_10 - 0.0000000000000001, mod_10 + gain_10); //C++ double should have a floating-point precision of up to 15 digits
    EXPECT_GE(mod_prime_10 + 0.0000000000000001, mod_10 + gain_10);
}
}