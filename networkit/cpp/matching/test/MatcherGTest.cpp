/*
 * MatcherGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/DibapGraphReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>
#include <networkit/matching/BMatcher.hpp>
#include <networkit/matching/BMatching.hpp>
#include <networkit/matching/BSuitorMatcher.hpp>
#include <networkit/matching/DynamicBSuitorMatcher.hpp>
#include <networkit/matching/LocalMaxMatcher.hpp>
#include <networkit/matching/Matcher.hpp>
#include <networkit/matching/Matching.hpp>
#include <networkit/matching/PathGrowingMatcher.hpp>
#include <networkit/matching/SuitorMatcher.hpp>
#include <networkit/viz/GraphLayoutAlgorithm.hpp>

namespace NetworKit {

class MatcherGTest : public testing::Test {
protected:
    Graph generateRandomWeightedGraph(count n) {
        // Generates a random undirected Graph with n nodes, a 50% chance for an edge between a
        // pair of nodes and without self-loops. Sets a random weight for each edge.
        auto G = ErdosRenyiGenerator(n, 0.5, false, false).generate();
        G = GraphTools::toWeighted(G);
        G.forEdges([&G](node u, node v) { G.setWeight(u, v, Aux::Random::integer(1, 99)); });
        return G;
    }

    // TODO use template instead of an overload when Matching base class is done
    bool hasUnmatchedNeighbors(const Graph &G, const BMatching &M) {
        for (const auto e : G.edgeRange())
            if (M.isUnmatched(e.u) && M.isUnmatched(e.v))
                return true;
        return false;
    }

    bool hasUnmatchedNeighbors(const Graph &G, const Matching &M) {
        for (const auto e : G.edgeRange())
            if (!M.isMatched(e.u) && !M.isMatched(e.v))
                return true;
        return false;
    }
};

TEST_F(MatcherGTest, testLocalMaxMatching) {
    {
        Graph G(10, true, true);
        EXPECT_THROW(LocalMaxMatcher{G}, std::runtime_error);
    }

    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v); });

    LocalMaxMatcher localMaxMatcher(G);

    TRACE("Start localMax matching");
    localMaxMatcher.run();
    Matching M = localMaxMatcher.getMatching();
    TRACE("Finished localMax matching");

    count numExpEdges = n / 2;
    bool isProper = M.isProper(G);
    EXPECT_TRUE(isProper);
    EXPECT_EQ(M.size(G), numExpEdges);

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
    DibapGraphReader reader;
    Graph airfoil1 = reader.read("input/airfoil1.gi");
    LocalMaxMatcher lmm(airfoil1);
    lmm.run();
    M = lmm.getMatching();
    isProper = M.isProper(airfoil1);
    EXPECT_TRUE(isProper);
    DEBUG("LocalMax on airfoil1 produces matching of size: ", M.size(G));
#endif
}

TEST_F(MatcherGTest, testLocalMaxMatchingDirectedWarning) {
    Graph G(2, false, true);
    G.addEdge(0, 1);
    EXPECT_THROW(LocalMaxMatcher localMaxMatcher(G), std::runtime_error);
}

TEST_F(MatcherGTest, testPgaMatchingOnWeightedGraph) {
    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v, Aux::Random::real()); });
    PathGrowingMatcher pgaMatcher(G);
    EXPECT_NO_THROW(pgaMatcher.run());
}

TEST_F(MatcherGTest, testPgaMatchingWithSelfLoops) {
    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v, Aux::Random::real()); });
    G.forNodes([&](node u) { G.addEdge(u, u); });
    EXPECT_THROW(PathGrowingMatcher pgaMatcher(G), std::invalid_argument);
}

TEST_F(MatcherGTest, testPgaMatching) {
    count n = 50;
    Graph G(n);
    G.forNodePairs([&](node u, node v) { G.addEdge(u, v); });
    PathGrowingMatcher pgaMatcher(G);

    DEBUG("Start PGA matching on 50-clique");

    pgaMatcher.run();
    Matching M = pgaMatcher.getMatching();

    count numExpEdges = n / 2;
    bool isProper = M.isProper(G);
    EXPECT_TRUE(isProper);
    EXPECT_EQ(M.size(G), numExpEdges);
    DEBUG("Finished PGA matching on 50-clique");

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
    DibapGraphReader reader;
    Graph airfoil1 = reader.read("input/airfoil1.gi");
    PathGrowingMatcher pga2(airfoil1);
    pga2.run();
    M = pga2.getMatching();
    isProper = M.isProper(airfoil1);
    EXPECT_TRUE(isProper);
    DEBUG("PGA on airfoil1 produces matching of size: ", M.size(G));
#endif
}

TEST_F(MatcherGTest, debugValidMatching) {
    METISGraphReader reader;
    Graph G = reader.read("coAuthorsDBLP.graph");

    LocalMaxMatcher pmatcher(G);
    pmatcher.run();
    Matching M = pmatcher.getMatching();

    bool isProper = M.isProper(G);
    EXPECT_TRUE(isProper);
}

TEST_F(MatcherGTest, testSuitorMatcher) {
    { // Directed graphs are not supported
        Graph G(10, true, true);
        EXPECT_THROW(SuitorMatcher{G}, std::runtime_error);
    }

    { // Graphs with self loops are not supported
        Graph G(10);
        G.addEdge(0, 0);
        G.addEdge(0, 0);
        EXPECT_THROW(SuitorMatcher{G}, std::runtime_error);
    }

    const edgeweight maxWeight = 10;

    const auto doTest = [&, maxWeight](Graph &G) -> void {
        // Test suitor matcher
        SuitorMatcher sm(G, false, true);
        sm.run();
        const auto M1 = sm.getMatching();
        EXPECT_TRUE(M1.isProper(G));
        EXPECT_FALSE(hasUnmatchedNeighbors(G, M1));

        GraphTools::sortEdgesByWeight(G, true);
        {
            auto G1 = G;
            G1.addEdge(0, 1, maxWeight);
            G1.addEdge(0, 2, maxWeight);
            EXPECT_THROW(SuitorMatcher(G1, true, true), std::runtime_error);
        }

        // Test sort suitor matcher
        SuitorMatcher ssm(G, true, true);
        ssm.run();
        const auto M2 = ssm.getMatching();
        EXPECT_TRUE(M2.isProper(G));
        EXPECT_FALSE(hasUnmatchedNeighbors(G, M2));

        // Matchings must be the same
        G.forNodes([&M1, &M2](node u) { EXPECT_EQ(M1.mate(u), M2.mate(u)); });
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);

        // Test unweighted
        auto G = METISGraphReader{}.read("input/PGPgiantcompo.graph");
        G.removeSelfLoops();
        G.removeMultiEdges();
        doTest(G);

        // Test weighted
        G = GraphTools::toWeighted(G);
        G.forEdges(
            [&G, maxWeight](node u, node v) { G.setWeight(u, v, Aux::Random::real(maxWeight)); });
        doTest(G);
    }
}

TEST_F(MatcherGTest, testBSuitorMatcherInvalidGraphDirected) {
    Graph G(10, true, true);
    EXPECT_THROW(BSuitorMatcher(G, 2), std::runtime_error);
}

TEST_F(MatcherGTest, testBSuitorMatcherInvalidGraphSelfLoops) {
    Graph G(10);
    G.addEdge(0, 0);
    G.addEdge(0, 0);
    EXPECT_THROW(BSuitorMatcher(G, 2), std::runtime_error);
}

TEST_F(MatcherGTest, testBSuitorMatcherTieBreaking) {
    auto G = METISGraphReader{}.read("input/tie.graph");
    G.removeSelfLoops();
    G.removeMultiEdges();

    BSuitorMatcher bsm(G, 1);
    bsm.run();
    bsm.buildBMatching();
    const auto M = bsm.getBMatching();

    EXPECT_TRUE(M.isProper(G));
    EXPECT_FALSE(hasUnmatchedNeighbors(G, M));
}

TEST_F(MatcherGTest, testBSuitorMatcherEqualsSuitorMatcher) {
    auto G = METISGraphReader{}.read("input/lesmis.graph");
    G.removeSelfLoops();
    G.removeMultiEdges();

    SuitorMatcher sm(G, false, false);
    sm.run();
    const auto M = sm.getMatching();

    BSuitorMatcher bsm(G, 1);
    bsm.run();
    bsm.buildBMatching();
    const auto bM = bsm.getBMatching();

    EXPECT_TRUE(bM.isProper(G));
    EXPECT_TRUE(M.isProper(G));
    EXPECT_FALSE(hasUnmatchedNeighbors(G, M));
    EXPECT_FALSE(hasUnmatchedNeighbors(G, bM));
}

TEST_F(MatcherGTest, testBSuitorMatcherConstantB) {
    for (int b : {2, 3, 4, 5}) {
        auto G = METISGraphReader{}.read("input/lesmis.graph");
        G.removeSelfLoops();
        G.removeMultiEdges();
        BSuitorMatcher bsm(G, b);
        bsm.run();
        bsm.buildBMatching();
        const auto M = bsm.getBMatching();
        EXPECT_TRUE(M.isProper(G));
        EXPECT_FALSE(hasUnmatchedNeighbors(G, M));
    }
}

TEST_F(MatcherGTest, testBSuitorMatcherDifferentB) {
    Aux::Random::setSeed(1, true);

    auto G = METISGraphReader{}.read("input/lesmis.graph");
    G.removeSelfLoops();
    G.removeMultiEdges();
    std::vector<count> b;
    for (count i = 0; i < G.numberOfNodes(); i++) {
        b.emplace_back(Aux::Random::integer(1, (G.numberOfNodes() - 1)));
    }

    BSuitorMatcher bsm(G, b);
    bsm.run();
    bsm.buildBMatching();
    const auto M = bsm.getBMatching();
    EXPECT_TRUE(M.isProper(G));
    EXPECT_FALSE(hasUnmatchedNeighbors(G, M));
}

TEST_F(MatcherGTest, testDynBSuitorInsertEdges) {
    for (int i = 0; i < 10000; i++) {
        // Aux::Random::setSeed(i, true);

        auto G = generateRandomWeightedGraph(100);
        std::vector<WeightedEdge> edges;
        std::vector<edgeweight> edgeWeights;
        count m = 10;
        // Select m edges of the graph, remove them but put them into edges for later insertion.
        // This will make sure that the graph is valid.
        for (auto j = 0; j < m; j++) {
            const auto [u, v] = GraphTools::randomEdge(G);
            assert(G.hasEdge(u, v));
            edges.emplace_back(u, v, G.weight(u, v));
            edgeWeights.emplace_back(G.weight(u, v));
            G.removeEdge(u, v);
        }

        const count b = 6;
        DynamicBSuitorMatcher dbsm(G, b);
        dbsm.run();

        for (auto &edge : edges) {
            G.addEdge(edge.u, edge.v, edge.weight);
        }

        dbsm.addEdges(edges);
        dbsm.buildBMatching();
        const auto dm = dbsm.getBMatching();
        auto dres = dm.getMatches();
        const auto dwm = dm.weight(G);

        BSuitorMatcher bsm(G, b);
        bsm.run();
        bsm.buildBMatching();
        const auto sm = bsm.getBMatching();
        auto res = sm.getMatches();
        const auto wm = sm.weight(G);

        if (dwm != wm) {

            for (size_t i = 0; i < edges.size(); i++) {
                G.removeEdge(edges[i].u, edges[i].v);
            }

            DynamicBSuitorMatcher dbsm2(G, b);
            dbsm2.run();

            for (size_t i = 0; i < edges.size(); i++) {
                G.addEdge(edges[i].u, edges[i].v, edgeWeights[i]);
            }

            Aux::Log::setLogLevel("INFO");

            dbsm2.addEdges(edges);
            dbsm2.buildBMatching();
            for (auto &myEdge : edges) {
                INFO("Edge: ", myEdge.u, ",", myEdge.v, ",", myEdge.weight);
            }

            for (size_t i = 0; i < res.size(); i++) {
                INFO("Node: ", i, " Identical: ", res.at(i) == dres.at(i), " ", res.at(i), " vs. ",
                     dres.at(i));
            }

            G.forNodes([&](node u) {
                std::string curN = std::to_string(u);
                std::string neighbors = "";
                G.forNeighborsOf(u, [&](node v) { neighbors += std::to_string(v) + ","; });
                INFO("Node ", curN, ": ", neighbors);
            });
            Aux::Log::setLogLevel("QUIET");

            NetworkitBinaryWriter nkWriter = NetworkitBinaryWriter{};
            nkWriter.write(G, "./G_inserted.nkbg");
        }

        // TODO: Change back to EXPECT in final code
        ASSERT_EQ(dwm, wm);
    }
}

TEST_F(MatcherGTest, testDynBSuitorRemoveEdges) {
    //
    for (int i = 0; i < 10000; i++) {
        // Aux::Random::setSeed(i, true);
        // Aux::Log::setLogLevel("INFO");

        auto G = generateRandomWeightedGraph(100);

        const count b = 6;
        DynamicBSuitorMatcher dbsm(G, b);
        dbsm.run();

        std::vector<Edge> edges;
        std::vector<edgeweight> edgeWeights;
        count m = 10;
        for (auto j = 0; j < m; j++) {
            const auto [u, v] = GraphTools::randomEdge(G);
            edges.emplace_back(u, v, G.weight(u, v));
            edgeWeights.emplace_back(G.weight(u, v));
            G.removeEdge(u, v);
        }

        dbsm.removeEdges(edges);
        dbsm.buildBMatching();
        const auto dm = dbsm.getBMatching();
        auto dres = dm.getMatches();
        const auto dwm = dm.weight(G);

        BSuitorMatcher bsm(G, b);
        bsm.run();
        bsm.buildBMatching();
        const auto sm = bsm.getBMatching();
        auto res = sm.getMatches();
        const auto wm = sm.weight(G);

        if (dwm != wm) {

            for (size_t i = 0; i < edges.size(); i++) {
                G.addEdge(edges[i].u, edges[i].v, edgeWeights[i]);
            }

            DynamicBSuitorMatcher dbsm2(G, b);
            dbsm2.run();

            for (auto &edge : edges) {
                G.removeEdge(edge.u, edge.v);
            }

            Aux::Log::setLogLevel("INFO");

            dbsm2.removeEdges(edges);
            dbsm2.buildBMatching();
            const auto dm2 = dbsm2.getBMatching();
            auto dres2 = dm2.getMatches();

            for (auto &myEdge : edges) {
                INFO("Edge: ", myEdge.u, ",", myEdge.v);
            }

            for (size_t i = 0; i < res.size(); i++) {
                if (res.at(i) != dres2.at(i)) {
                    INFO("Node: ", i, " Identical: ", res.at(i) == dres2.at(i), " ", res.at(i),
                         " vs. ", dres2.at(i));
                }
            }

            G.forNodes([&](node u) {
                std::string curN = std::to_string(u);
                std::string neighbors = "";
                G.forNeighborsOf(u, [&](node v) { neighbors += std::to_string(v) + ","; });
                INFO("Node ", curN, ": ", neighbors);
            });
            Aux::Log::setLogLevel("QUIET");
        }

        // TODO: change back to expect in final code
        ASSERT_EQ(dwm, wm);
    }
}
TEST_F(MatcherGTest, testDynBSuitorMulti) {
    Aux::Log::setLogLevel("INFO");
    NetworKit::Graph G = NetworKit::Graph(10, true, false);

    G.addEdge(0, 4, 1);
    G.addEdge(0, 5, 1);
    G.addEdge(0, 6, 1);
    G.addEdge(0, 7, 21);
    G.addEdge(0, 8, 30);

    G.addEdge(1, 2, 1);
    G.addEdge(1, 4, 1);
    G.addEdge(1, 6, 1);
    G.addEdge(1, 7, 56);
    G.addEdge(1, 9, 93);

    G.addEdge(2, 4, 59);
    G.addEdge(2, 5, 70);
    G.addEdge(2, 8, 1);
    G.addEdge(2, 9, 65);

    G.addEdge(3, 4, 1);
    G.addEdge(3, 6, 1);
    G.addEdge(3, 8, 50);

    G.addEdge(4, 6, 60);
    G.addEdge(4, 7, 48);

    G.addEdge(7, 8, 40);

    // New test edges for finding loops
    // G.addEdge(4, 6, 40);
    // G.addEdge(2, 6, 1);
    // G.addEdge(3, 6, 2);

    BSuitorMatcher bsm(G, 2);
    bsm.run();
    bsm.buildBMatching();
    auto sm = bsm.getBMatching();

    auto res = sm.getMatches();

    // auto idx = 0;
    // INFO("Initial matching:");
    // for (auto &u : res) {
    //     INFO(idx, ": ");
    //     for (auto v : u) {
    //         INFO(v, ", ");
    //     }
    //     idx++;
    // }

    NetworKit::DynamicBSuitorMatcher dynBMatcher(G, 2);
    dynBMatcher.run();
    dynBMatcher.buildBMatching();

    auto bMatching = dynBMatcher.getBMatching();
    auto bres = bMatching.getMatches();

    for (size_t i = 0; i < res.size(); i++) {
        INFO(res.at(i), " vs. ", bres.at(i));
    }

    G.addEdge(7, 9, 84);

    BSuitorMatcher bsm2(G, 2);
    bsm2.run();
    bsm2.buildBMatching();
    auto sm2 = bsm2.getBMatching();

    auto res2 = sm2.getMatches();

    NetworKit::WeightedEdge newEdge(7, 9, 84);
    dynBMatcher.addEdge(newEdge);
    dynBMatcher.buildBMatching();
    auto bMatching2 = dynBMatcher.getBMatching();
    auto bres2 = bMatching2.getMatches();

    for (size_t i = 0; i < res2.size(); i++) {
        INFO(res2.at(i), " vs. ", bres2.at(i));
    }

    G.removeEdge(7, 9);
    dynBMatcher.removeEdge(newEdge);
    dynBMatcher.buildBMatching();
    auto bMatching3 = dynBMatcher.getBMatching();
    auto bres3 = bMatching3.getMatches();

    for (size_t i = 0; i < res.size(); i++) {
        INFO(res.at(i), " vs. ", bres3.at(i));
    }

    // for (size_t i = 0; i < dynBMatcher.Proposed.size(); i++) {
    //     INFO(*bsm2.Proposed.at(i), " vs. ", *dynBMatcher.Proposed.at(i));
    // }
}

TEST_F(MatcherGTest, testDynBSuitorUnsaturated) {
    Aux::Log::setLogLevel("INFO");
    NetworKit::Graph G = NetworKit::Graph(6, true, false);

    G.addEdge(0, 1, 1);
    G.addEdge(0, 2, 1);

    G.addEdge(1, 2, 1);
    G.addEdge(1, 3, 1);

    G.addEdge(3, 4, 9);

    G.addEdge(4, 5, 6);

    BSuitorMatcher bsm(G, 2);
    bsm.run();
    bsm.buildBMatching();
    auto sm = bsm.getBMatching();

    auto res = sm.getMatches();

    // auto idx = 0;
    // INFO("Initial matching:");
    // for (auto &u : res) {
    //     INFO(idx, ": ");
    //     for (auto v : u) {
    //         INFO(v, ", ");
    //     }
    //     idx++;
    // }

    NetworKit::DynamicBSuitorMatcher dynBMatcher(G, 2);
    dynBMatcher.run();
    dynBMatcher.buildBMatching();

    auto bMatching = dynBMatcher.getBMatching();
    auto bres = bMatching.getMatches();

    for (size_t i = 0; i < res.size(); i++) {
        INFO(res.at(i), " vs. ", bres.at(i));
    }

    G.addEdge(3, 5, 4);

    BSuitorMatcher bsm2(G, 2);
    bsm2.run();
    bsm2.buildBMatching();
    auto sm2 = bsm2.getBMatching();

    auto res2 = sm2.getMatches();

    NetworKit::WeightedEdge newEdge(3, 5, 4);
    dynBMatcher.addEdge(newEdge);
    dynBMatcher.buildBMatching();
    auto bMatching2 = dynBMatcher.getBMatching();
    auto bres2 = bMatching2.getMatches();

    for (size_t i = 0; i < res2.size(); i++) {
        INFO(res2.at(i), " vs. ", bres2.at(i));
    }

    // for (size_t i = 0; i < dynBMatcher.Proposed.size(); i++) {
    //     INFO(*bsm2.Proposed.at(i), " vs. ", *dynBMatcher.Proposed.at(i));
    // }
}

} // namespace NetworKit
