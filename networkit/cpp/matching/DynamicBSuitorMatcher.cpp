#include <cassert>
#include <chrono>

#include <networkit/matching/DynamicBSuitorMatcher.hpp>

namespace NetworKit {

void DynamicBSuitorMatcher::processEdgeInsertion(const WeightedEdge &edge) {

    G->forNodes([&](node u) { Suitors.at(u)->reset(); });

    // affected[edge.u] = true;
    affectedNodes.clear();
    findAffectedNodes(edge.u, edge.v, Operation::Insert);
    size_t batchId = batchTracker[std::make_pair(edge.u,edge.v)];
    updateAffectedNodes(batchId);
    // affected[edge.u] = false;

    // affected[edge.v] = true;
    // affectedNodes.clear();
    // findAffectedNodes(edge.v, edge.u, Operation::Insert);
    // updateAffectedNodes(batchId);

    // affected[edge.u] = affected[edge.v] = false;
}

void DynamicBSuitorMatcher::processEdgeInsertionNew(const WeightedEdge &edge) {

    node u = edge.u;
    node v = edge.v;
    // INFO("Start processing ", u, ",", v);
    edgeweight w = G->weight(u,v);
    
    auto edgeHash = u < v ? std::make_pair(u, v) : std::make_pair(v, u);
    size_t batchId = batchTracker[edgeHash];

    DynBNode startU = Suitors.at(u)->insert({v,w});
    DynBNode startV = Suitors.at(v)->insert({u,w});

    // INFO("StartU: ", startU.id);
    // INFO("StartV: ", startV.id);

    if (startU.id != none) {
        Suitors.at(startU.id)->remove(u);
    }

    if (startV.id != none) {
        Suitors.at(startV.id)->remove(v);
    }

    if(startU.id != none) {
        trackUpdatePath(batchId, startU.id);
    }

    if(startV.id != none) {
        trackUpdatePath(batchId, startV.id);
    }
}

void DynamicBSuitorMatcher::trackUpdatePath(size_t batchId, node start, bool recursiveCall) {
    // INFO("Start node: ", start);
    std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());

    bool done = false;

    node current = start;
    node partner = Suitors.at(current)->min.id;
    auto heaviest = Suitors.at(current)->min.weight;
    // INFO("Starting findAffectedNodes8 for id ", batchId, " with: \n cur: ", current, ", partner: ", partner,
         ", weight: ", heaviest);

    edgeweight prev = std::numeric_limits<edgeweight>::max();

    std::vector<DynBNode> looseEnds;
    // int numIter = 0;

    do {
        done = true;

        G->forNeighborsOf(current, [&](node x, edgeweight weight) {
            auto edgeHash = current < x ? std::make_pair(current, x) : std::make_pair(x, current);
            // ignore edges that that still need to be processed
            if (batchTracker.contains(edgeHash)) {
                // INFO("Processing batch item: ", batchId, " Edge ", current, ",", x, " has id ", batchTracker[edgeHash]);
                if (batchId < batchTracker[edgeHash])
                    return;
            }

            if (Suitors.at(current)->hasPartner(x))
                return;

            const auto z = Suitors.at(x)->min;
            // INFO("Node: ", current, " (weight:", heaviest, ") Evaluating ", x, " (edge_weight: ", weight, " min: ", z.id,
                 ", min_weight: ", z.weight, ")");

            if ((weight > heaviest || (weight == heaviest && x < partner))
                && (weight > z.weight || (weight == z.weight && current < z.id))
                && (weight <= prev)) {
                partner = x;
                heaviest = weight;
                return;
            }
        });

        // INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
            //  ", weight: ", heaviest);

        // line 10-12
        if (partner == none || heaviest < Suitors.at(partner)->min.weight
            || (heaviest == Suitors.at(partner)->min.weight
                && current > Suitors.at(partner)->min.id)) {
            break;
        }

        // INFO("Suitors at current ", current, ": ", *Suitors.at(current));
        // INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
        // INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
        //      " (weight: ", (Suitors.at(partner)->min).weight, ")");

        // INFO("Inserting suitor ", current, " into ", partner, " and vice versa.");
        DynBNode prevCurrent = Suitors.at(current)->insert({partner, heaviest});
        DynBNode prevPartner = Suitors.at(partner)->insert({current, heaviest});
        affectedNodesPerRun++;

        if(prevCurrent.id != none) {
            // INFO("Current was saturated. Removing current from prevCurrent ", prevCurrent.id);
            Suitors.at(prevCurrent.id)->remove(current);
            looseEnds.emplace_back(DynBNode{prevCurrent.id,Suitors.at(prevCurrent.id)->min.weight});
        }


        if (prevPartner.id != none) {
            // INFO("Removing suitor ", prevPartner.id, " from ", partner);
            // INFO("Suitors at y ", prevPartner.id, ": ", *Suitors.at(prevPartner.id));

            // INFO("Removing suitor ", partner, " from ", prevPartner.id);
            Suitors.at(prevPartner.id)->remove(partner);
            // affected[y.id] = true;
            current = prevPartner.id;
            done = false;
        }

        // line 21-23
        prev = heaviest;
        partner = Suitors.at(current)->min.id;
        heaviest = Suitors.at(current)->min.weight;
        // numIter++;
    } while (!done);

    // if (looseEnds.size() > 0) {
    //     INFO("Recursive: ",recursiveCall ," Timestamp: " , ms.count(), " -> start: " ,start , " iter: " ,numIter , " loose ends: " , looseEnds.size() ," (" , looseEnds[0].id , ")" );
    // }
    for (auto &looseEnd : looseEnds) {
        
        trackUpdatePath(batchId, looseEnd.id, true);
    }
}

void DynamicBSuitorMatcher::processEdgeRemovalNew(const Edge &edge) {

    node u = edge.u;
    node v = edge.v;
    // INFO("Start processing ", u, ",", v);
    
    auto edgeHash = u < v ? std::make_pair(u, v) : std::make_pair(v, u);
    size_t batchId = batchTracker[edgeHash];

    Suitors.at(u)->remove(v);
    Suitors.at(v)->remove(u);

    trackUpdatePath(batchId, u);    
    trackUpdatePath(batchId, v);
}


void DynamicBSuitorMatcher::processEdgeRemoval(const Edge &edge) {

    G->forNodes([&](node u) { Suitors.at(u)->reset(); });

    Suitors.at(edge.u)->remove(edge.v);
    Suitors.at(edge.v)->remove(edge.u);

    // affected[edge.u] = affected[edge.v] = true;
    affectedNodes.clear();
    affectedNodes.emplace_back(edge.u);
    findAffectedNodes(edge.u, edge.v, Operation::Remove);
    size_t batchId = batchTracker[std::make_pair(edge.u,edge.v)];
    updateAffectedNodes(batchId);
    // affected[edge.u] = false;

    // affected[edge.v] = true;
    affectedNodes.clear();
    affectedNodes.emplace_back(edge.v);
    findAffectedNodes(edge.v, edge.u, Operation::Remove);
    updateAffectedNodes(batchId);
    // affected[edge.v] = false;
}

void DynamicBSuitorMatcher::addEdges(std::vector<WeightedEdge> &edges, bool sort) {

    affectedNodesPerRun = 0;
    batchTracker.clear();

    // auto convertToEdge = [](const WeightedEdge &weightedEdge) {
    //     return Edge(weightedEdge.u, weightedEdge.v);
    // };
    // std::transform(edges.begin(), edges.end(), std::back_inserter(edgeBatch), convertToEdge);
    if (sort) {
        std::sort(edges.begin(), edges.end(),
                [](const WeightedEdge &a, const WeightedEdge &b) { return a.weight > b.weight; });
    }

    for (size_t i = 0; i < edges.size(); i++)
    {
        auto edgeHash = edges[i].u < edges[i].v ? std::make_pair(edges[i].u, edges[i].v) : std::make_pair(edges[i].v, edges[i].u);
        batchTracker[edgeHash] = i;
        // INFO("Batch: ", batchTracker[edgeHash], " -> ", edges[i].u, ",", edges[i].v);
    }

    // for (std::unordered_map<std::pair<int,int>,uint8_t,PairHash>::const_iterator it = batchTracker.begin();
    //     it != batchTracker.end(); ++it) {
    //     INFO(" [" << it->first.first << "," << it->first.second << " -> " << it->second << "]";
    //     std::cout << std::endl;

    //     }


    for (const auto &edge : edges) {
        if ((Suitors.at(edge.u)->hasPartner(edge.v) && Suitors.at(edge.v)->hasPartner(edge.u))
            || !isBetterMatch(edge.u, edge.v, edge.weight)
            || !isBetterMatch(edge.v, edge.u, edge.weight)) {
            // INFO("Edge ", edge.u, ",", edge.v,
            //      " ignored, since min(u): ", Suitors.at(edge.u)->min.weight,
            //      " min(v): ", Suitors.at(edge.v)->min.weight);
            continue;
        }
        // INFO("Edge ", edge.u, ",", edge.v, " better choice. Min at ", edge.u, ": ", Suitors.at(edge.u)->min.weight, " Min at ", edge.v, ": ", Suitors.at(edge.v)->min.weight);
        // affectedNodes.clear();

        // processEdgeInsertion(edge);
        processEdgeInsertionNew(edge);
        // affected[edge.u] = affected[edge.v] = false;

        // #ifndef NDEBUG
        //         G->forNodes([&](node u) {
        //             for (auto s : Suitors.at(u)->partners) {
        //                 assert(Suitors.at(s.id)->hasPartner(u));
        //             }
        //         });
        // #endif
    }
}

void DynamicBSuitorMatcher::removeEdges(std::vector<Edge> &edges) {

    affectedNodesPerRun = 0;
    batchTracker.clear();

    for (size_t i = 0; i < edges.size(); i++)
    {
        auto edgePair = std::make_pair(edges[i].u, edges[i].v);
        batchTracker[edgePair] = i;
    }


    // edgeBatch = edges;

    for (const auto &edge : edges) {
        assert(!G->hasEdge(edge.u, edge.v));
        if (Suitors.at(edge.u)->hasPartner(edge.v)) {

            // affectedNodes.clear();
            // assert(numberOfAffectedEquals(0));

            // processEdgeRemoval(edge);
            processEdgeRemovalNew(edge);
            // affected[edge.u] = affected[edge.v] = false;

            // #ifndef NDEBUG
            //             G->forNodes([&](node u) {
            //                 for (auto s : Suitors.at(u)->partners) {
            //                     assert(Suitors.at(s.id)->hasPartner(u));
            //                 }
            //             });
            // #endif
        }
    }
}

// void DynamicBSuitorMatcher::findAffectedNodes(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//             // ignore edges that that still need to be processed
//             if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) !=
//             edgeBatch.end()
//                 || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                        != edgeBatch.end())
//                 return;
//             if (Suitors.at(current)->hasPartner(x))
//                 return;

//             const auto z = Suitors.at(x)->min;
//             if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
//                 if ((weight > heaviest || (weight == heaviest && x < partner))
//                     && (weight > z.weight || (weight == z.weight && current < z.id))
//                     && (weight <= prev)) {
//                     partner = x;
//                     heaviest = weight;
//                     return;
//                 }
//             }
//         });
//         done = true;
//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         if (partner == none) {
//             affected[current] = false;
//             done = true;
//             continue;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         const auto ps = Suitors.at(partner)->min;
//         if (heaviest > ps.weight || (heaviest == ps.weight && current < ps.id)) {
//             const auto y = Suitors.at(partner)->popMinIfFull();

//             if (y.id != none) {
//                 INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//                 INFO("Removing suitor ", y.id, " from ", partner);
//             }

//             INFO("Inserting suitor ", current, " into ", partner);
//             DynBNode newS = Suitors.at(partner)->insert({current, heaviest});
//             // if (newS.id != none) {
//             //     Suitors.at(newS.id)->remove(partner);
//             // }
//             affected[partner] = true;
//             affectedNodes.emplace_back(partner);

//             if (y.id != none) {
//                 INFO("Removing suitor ", partner, " from ", y.id);
//                 Suitors.at(y.id)->remove(partner);
//                 affected[y.id] = true;
//                 current = y.id;
//                 done = false;
//             }
//         } else {
//             affected[current] = false;
//         }
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//         INFO("End of current loop: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);
//     } while (!done);
// }

// void DynamicBSuitorMatcher::findAffectedNodes2(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes2 with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         // line 8
//         done = true;

//         // line 9
//         // if (op == Operation::Insert) {
//         G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//             // ignore edges that that still need to be processed
//             if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) !=
//             edgeBatch.end()
//                 || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                        != edgeBatch.end())
//                 return;
//             if (Suitors.at(current)->hasPartner(x))
//                 return;

//             const auto z = Suitors.at(x)->min;
//             INFO("Node: ", current, "(weight:", heaviest, ") Evaluating ", x, " (min: ",
//             z.weight,
//                  ", weight: ", weight, ")");
//             // !affected[x] &&
//             if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
//                 if ((weight > heaviest || (weight == heaviest && x < partner))
//                     && (weight > z.weight || (weight == z.weight && current < z.id))
//                     && (weight <= prev)) {
//                     partner = x;
//                     heaviest = weight;
//                     return;
//                 }
//             }
//         });

//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         // line 10-12
//         if (partner == none || heaviest < Suitors.at(partner)->min.weight
//             || (heaviest == Suitors.at(partner)->min.weight
//                 && current > Suitors.at(partner)->min.id)) {
//             affected[current] = false;
//             break;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         // line 13-15
//         INFO("Inserting suitor ", current, " into ", partner);
//         DynBNode y = Suitors.at(partner)->insert({current, heaviest});
//         affected[partner] = true;
//         affectedNodes.emplace_back(partner);

//         if (y.id != none) {
//             INFO("Removing suitor ", y.id, " from ", partner);
//         }

//         // line 16-20
//         if (y.id != none) {
//             INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

//             INFO("Removing suitor ", partner, " from ", y.id);
//             Suitors.at(y.id)->remove(partner);
//             affected[y.id] = true;
//             current = y.id;
//             done = false;
//         }

//         // line 21-23
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//     } while (!done);
// }

// void DynamicBSuitorMatcher::findAffectedNodes3(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes3 with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         // line 8
//         done = true;

//         // line 9
//         G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//             // ignore edges that that still need to be processed
//             if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) !=
//             edgeBatch.end()
//                 || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                        != edgeBatch.end())
//                 return;
//             if (Proposed.at(current)->hasPartner(x))
//                 return;

//             const auto z = Suitors.at(x)->min;
//             INFO("Node: ", current, "(weight:", heaviest, ") Evaluating ", x, " (min: ",
//             z.weight,
//                  ", weight: ", weight, ")");
//             // !affected[x] &&
//             if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
//                 if ((weight > heaviest || (weight == heaviest && x < partner))
//                     && (weight > z.weight || (weight == z.weight && current < z.id))
//                     && (weight <= prev)) {
//                     partner = x;
//                     heaviest = weight;
//                     return;
//                 }
//             }
//         });

//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         // line 10-12
//         if (partner == none || heaviest < Suitors.at(partner)->min.weight
//             || (heaviest == Suitors.at(partner)->min.weight
//                 && current > Suitors.at(partner)->min.id)) {
//             affected[current] = false;
//             break;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         // line 13-15
//         INFO("Inserting suitor ", current, " into ", partner);
//         // T-invariant
//         Proposed.at(current)->insert({partner, heaviest});
//         DynBNode y = Suitors.at(partner)->insert({current, heaviest});
//         affected[partner] = true;
//         affectedNodes.emplace_back(partner);

//         if (y.id != none) {
//             INFO("Removing suitor ", y.id, " from ", partner);
//         }

//         // line 16-20
//         if (y.id != none) {
//             INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

//             INFO("Removing suitor ", partner, " from ", y.id);
//             // T-invariant
//             Proposed.at(partner)->remove(y.id);
//             Proposed.at(y.id)->remove(partner);
//             Suitors.at(y.id)->remove(partner);
//             affected[y.id] = true;
//             current = y.id;
//             done = false;
//         }

//         // line 21-23
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//     } while (!done);
// }

// void DynamicBSuitorMatcher::findAffectedNodes4(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes4 with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         // line 8
//         done = true;

//         // line 9
//         G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//             // ignore edges that that still need to be processed
//             if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) !=
//             edgeBatch.end()
//                 || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                        != edgeBatch.end())
//                 return;
//             if (Proposed.at(current)->hasPartner(x))
//                 return;

//             const auto z = Suitors.at(x)->min;
//             INFO("Node: ", current, "(weight:", heaviest, ") Evaluating ", x, " (min: ",
//             z.weight,
//                  ", weight: ", weight, ")");
//             // if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
//             if ((weight > heaviest || (weight == heaviest && x < partner))
//                 && (weight > z.weight || (weight == z.weight && current < z.id))
//                 && (weight <= prev)) {
//                 partner = x;
//                 heaviest = weight;
//                 return;
//             }
//             // }
//         });

//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         // line 10-12
//         if (partner == none || heaviest < Suitors.at(partner)->min.weight
//             || (heaviest == Suitors.at(partner)->min.weight
//                 && current > Suitors.at(partner)->min.id)) {
//             affected[current] = false;
//             break;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         // line 13-15
//         INFO("Inserting suitor ", current, " into ", partner);
//         // T-invariant
//         Proposed.at(current)->insert({partner, heaviest});
//         DynBNode y = Suitors.at(partner)->insert({current, heaviest});
//         Suitors.at(partner)->num_visits++;
//         affected[partner] = true;
//         affectedNodes.emplace_back(partner);

//         if (y.id != none) {
//             INFO("Removing suitor ", y.id, " from ", partner);
//         }

//         // line 16-20
//         if (y.id != none) {
//             INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

//             INFO("Removing suitor ", partner, " from ", y.id);
//             // T-invariant
//             Proposed.at(partner)->remove(y.id);
//             Proposed.at(y.id)->remove(partner);
//             Suitors.at(y.id)->remove(partner);
//             Suitors.at(y.id)->num_visits++;
//             affected[y.id] = true;
//             current = y.id;
//             done = false;
//         }

//         // line 21-23
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//     } while (!done);
// }

// void DynamicBSuitorMatcher::findAffectedNodes5(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes5 with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         // line 8
//         done = true;

//         // line 9
//         // if (op == Operation::Insert) {
//         bool skipNeighborCheck = false;
//         if (Suitors.at(current)->activeLooseEnd != none) {
//             skipNeighborCheck = true;
//             partner = Suitors.at(current)->activeLooseEnd;
//             heaviest = G->weight(current, partner);
//             Suitors.at(current)->activeLooseEnd = none;
//             if (heaviest < Suitors.at(partner)->min.weight
//                 || (heaviest == Suitors.at(partner)->min.weight
//                     && current > Suitors.at(partner)->min.id)) {
//                 skipNeighborCheck = false;
//                 partner = none;
//                 heaviest = 0;
//             }
//         }
//         if (!skipNeighborCheck) {
//             G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//                 // ignore edges that that still need to be processed
//                 if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x))
//                         != edgeBatch.end()
//                     || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                            != edgeBatch.end())
//                     return;
//                 if (Suitors.at(current)->hasPartner(x))
//                     return;

//                 const auto z = Suitors.at(x)->min;
//                 INFO("Node: ", current, "(weight:", heaviest, ") Evaluating ", x,
//                      " (min: ", z.weight, ", weight: ", weight, ")");
//                 // !affected[x] &&
//                 if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
//                     if ((weight > heaviest || (weight == heaviest && x < partner))
//                         && (weight > z.weight || (weight == z.weight && current < z.id))
//                         && (weight <= prev)) {
//                         partner = x;
//                         heaviest = weight;
//                         return;
//                     }
//                 }
//             });
//         }

//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         // line 10-12
//         if (partner == none || heaviest < Suitors.at(partner)->min.weight
//             || (heaviest == Suitors.at(partner)->min.weight
//                 && current > Suitors.at(partner)->min.id)) {
//             auto potentialLooseEnd = Suitors.at(current)->looseEnd;
//             if (potentialLooseEnd != none) {
//                 INFO("Adding loose end ", current, " to ", potentialLooseEnd);
//                 Suitors.at(potentialLooseEnd)->activeLooseEnd = current;
//                 Suitors.at(current)->looseEnd = none;
//             }

//             affected[current] = false;
//             break;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         // line 13-15
//         INFO("Inserting suitor ", current, " into ", partner);
//         DynBNode y = Suitors.at(partner)->insert({current, heaviest});
//         if (Suitors.at(partner)->activeLooseEnd == current) {
//             Suitors.at(partner)->activeLooseEnd = none;
//         }
//         affected[partner] = true;
//         affectedNodes.emplace_back(partner);

//         if (y.id != none) {
//             INFO("Removing suitor ", y.id, " from ", partner);
//         }

//         // line 16-20
//         if (y.id != none) {
//             INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

//             INFO("Removing suitor ", partner, " from ", y.id);
//             Suitors.at(y.id)->remove(partner);
//             affected[y.id] = true;
//             current = y.id;
//             done = false;
//         }

//         // line 21-23
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//     } while (!done);
// }

// void DynamicBSuitorMatcher::findAffectedNodes6(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes6 with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         // line 8
//         done = true;

//         // line 9
//         // if (op == Operation::Insert) {
//         G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//             // ignore edges that that still need to be processed
//             if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) !=
//             edgeBatch.end()
//                 || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                        != edgeBatch.end())
//                 return;
//             if (Suitors.at(current)->hasPartner(x))
//                 return;

//             const auto z = Suitors.at(x)->min;
//             INFO("Node: ", current, "(weight:", heaviest, ") Evaluating ", x, " (min: ",
//             z.weight,
//                  ", weight: ", weight, ")");
//             // !affected[x] &&
//             if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
//                 if ((weight > heaviest || (weight == heaviest && x < partner))
//                     && (weight > z.weight || (weight == z.weight && current < z.id))
//                     && (weight <= prev)) {
//                     partner = x;
//                     heaviest = weight;
//                     return;
//                 }
//             }
//         });

//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         // line 10-12
//         if (partner == none || heaviest < Suitors.at(partner)->min.weight
//             || (heaviest == Suitors.at(partner)->min.weight
//                 && current > Suitors.at(partner)->min.id)) {
//             affected[current] = false;
//             break;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         // line 13-15
//         INFO("Inserting suitor ", current, " into ", partner);
//         DynBNode y = Suitors.at(partner)->insert({current, heaviest});
//         affected[partner] = true;
//         affectedNodes.emplace_back(partner);

//         if (y.id != none) {
//             INFO("Removing suitor ", y.id, " from ", partner);
//         }

//         // line 16-20
//         if (y.id != none) {
//             INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

//             INFO("Removing suitor ", partner, " from ", y.id);
//             Suitors.at(y.id)->remove(partner);
//             affected[y.id] = true;
//             current = y.id;
//             done = false;
//         }

//         // line 21-23
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//     } while (!done);
// }

// void DynamicBSuitorMatcher::findAffectedNodes7(node u, node v, Operation op) {
//     bool done = false;

//     node current = u;
//     node partner = (op == Operation::Insert) ? v : none;
//     auto heaviest = G->weight(current, v);
//     INFO("Starting findAffectedNodes7 with: \n cur: ", current, ", partner: ", partner,
//          ", weight: ", heaviest);

//     edgeweight prev = std::numeric_limits<edgeweight>::max();

//     do {
//         // line 8
//         done = true;

//         // line 9
//         // if (op == Operation::Insert) {
//         G->forNeighborsOf(current, [&](node x, edgeweight weight) {
//             // ignore edges that that still need to be processed
//             if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) !=
//             edgeBatch.end()
//                 || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
//                        != edgeBatch.end())
//                 return;
//             if (Suitors.at(current)->hasPartner(x))
//                 return;

//             const auto z = Suitors.at(x)->min;
//             INFO("Node: ", current, "(weight:", heaviest, ") Evaluating ", x, " (min: ",
//             z.weight,
//                  ", weight: ", weight, ")");
//             // !affected[x] &&
//             if (!affected[x] || (op == Operation::Insert && (x == u || x == v))
//                 || (op == Operation::Remove && (x == u))) {
//                 if ((weight > heaviest || (weight == heaviest && x < partner))
//                     && (weight > z.weight || (weight == z.weight && current < z.id))
//                     && (weight <= prev)) {
//                     partner = x;
//                     heaviest = weight;
//                     return;
//                 }
//             }
//         });

//         INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
//              ", weight: ", heaviest);

//         // line 10-12
//         if (partner == none || heaviest < Suitors.at(partner)->min.weight
//             || (heaviest == Suitors.at(partner)->min.weight
//                 && current > Suitors.at(partner)->min.id)) {
//             affected[current] = false;
//             break;
//         }

//         INFO("Suitors at current ", current, ": ", *Suitors.at(current));
//         INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
//         INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
//              " (weight: ", (Suitors.at(partner)->min).weight, ")");

//         // line 13-15
//         INFO("Inserting suitor ", current, " into ", partner);
//         DynBNode y = Suitors.at(partner)->insert({current, heaviest});
//         affected[partner] = true;
//         affectedNodes.emplace_back(partner);

//         if (y.id != none) {
//             INFO("Removing suitor ", y.id, " from ", partner);
//         }

//         // line 16-20
//         if (y.id != none) {
//             INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

//             INFO("Removing suitor ", partner, " from ", y.id);
//             Suitors.at(y.id)->remove(partner);
//             affected[y.id] = true;
//             current = y.id;
//             done = false;
//         }

//         // line 21-23
//         prev = heaviest;
//         partner = Suitors.at(current)->min.id;
//         heaviest = Suitors.at(current)->min.weight;
//     } while (!done);
// }

void DynamicBSuitorMatcher::findAffectedNodes(node u, node v, Operation op) {
    bool done = false;

    node current = u;
    node partner = (op == Operation::Insert) ? v : Suitors.at(u)->min.id;
    auto heaviest = G->weight(current, partner);
    auto batchId = batchTracker[std::make_pair(current, partner)];
    INFO("Starting findAffectedNodes8 for id ", batchId, " with: \n cur: ", current, ", partner: ", partner,
         ", weight: ", heaviest);

    edgeweight prev = std::numeric_limits<edgeweight>::max();



    do {
        // line 8
        done = true;

        // line 9
        // if (op == Operation::Insert) {
        G->forNeighborsOf(current, [&](node x, edgeweight weight) {
            auto edgeHash = current < x ? std::make_pair(current, x) : std::make_pair(x, current);
            // ignore edges that that still need to be processed
            if (batchTracker.contains(edgeHash)) {
                INFO("Processing batch item: ", batchId, " Edge ", current, ",", x, " has id ", batchTracker[edgeHash]);
                if (batchId < batchTracker[edgeHash])
                    return;
            }
            // if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) != edgeBatch.end()
            //     || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
            //            != edgeBatch.end())
            //     return;
            if (Suitors.at(current)->hasPartner(x))
                return;

            const auto z = Suitors.at(x)->min;
            INFO("Node: ", current, " (weight:", heaviest, ") Evaluating ", x, " (edge_weight: ", weight, " min: ", z.id,
                 ", min_weight: ", z.weight, ")");
            // if (op == Operation::Remove && (z.id == u || z.id == v))
            //     return;
            // && (x == u || x == v)
            // if (!affected[x] || (op == Operation::Insert )
            //     || (op == Operation::Remove)) {
            if ((weight > heaviest || (weight == heaviest && x < partner))
                && (weight > z.weight || (weight == z.weight && current < z.id))
                && (weight <= prev)) {
                partner = x;
                heaviest = weight;
                return;
            }
            // }
        });

        INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
             ", weight: ", heaviest);

        // line 10-12
        if (partner == none || heaviest < Suitors.at(partner)->min.weight
            || (heaviest == Suitors.at(partner)->min.weight
                && current > Suitors.at(partner)->min.id)) {
            // affected[current] = false;
            break;
        }

        INFO("Suitors at current ", current, ": ", *Suitors.at(current));
        INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
        INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
             " (weight: ", (Suitors.at(partner)->min).weight, ")");

        // line 13-15
        INFO("Inserting suitor ", current, " into ", partner);
        DynBNode y = Suitors.at(partner)->insert({current, heaviest});
        // affected[partner] = true;
        affectedNodes.emplace_back(partner);

        if (y.id != none) {
            INFO("Removing suitor ", y.id, " from ", partner);
        }

        // line 16-20
        if (y.id != none) {
            INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));

            INFO("Removing suitor ", partner, " from ", y.id);
            Suitors.at(y.id)->remove(partner);
            // affected[y.id] = true;
            current = y.id;
            done = false;
        }

        // line 21-23
        prev = heaviest;
        partner = Suitors.at(current)->min.id;
        heaviest = Suitors.at(current)->min.weight;
    } while (!done);
}

// void DynamicBSuitorMatcher::updateAffectedNodes() {
//     INFO("Starting updateAffectedNodes");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         Suitors.at(y.id)->insert({x, y.weight});
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

// void DynamicBSuitorMatcher::updateAffectedNodes2() {
//     affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
//     INFO("Starting updateAffectedNodes2");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//         INFO("Inserting suitor ", x, " into ", y.id);
//         auto weight = G->weight(y.id, x);
//         auto minWeight = Suitors.at(y.id)->min.weight;
//         if (minWeight < weight || (minWeight == weight && Suitors.at(y.id)->min.id > x)) {
//             DynBNode newS = Suitors.at(y.id)->insert({x, y.weight});
//             if (newS.id != none) {
//                 Suitors.at(newS.id)->remove(y.id);
//             }
//         }
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

// void DynamicBSuitorMatcher::updateAffectedNodes3() {
//     affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
//     INFO("Starting updateAffectedNodes3");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//         INFO("Inserting suitor ", x, " into ", y.id);
//         auto weight = G->weight(y.id, x);
//         auto minWeight = Suitors.at(y.id)->min.weight;
//         if (minWeight < weight || (minWeight == weight && Suitors.at(y.id)->min.id > x)) {
//             Proposed.at(x)->insert({y.id, y.weight});
//             DynBNode newS = Suitors.at(y.id)->insert({x, y.weight});

//             if (newS.id != none) {
//                 Suitors.at(newS.id)->remove(y.id);
//                 Proposed.at(y.id)->remove(newS.id);
//             }
//         }
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

// void DynamicBSuitorMatcher::updateAffectedNodes4() {
//     affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
//     INFO("Starting updateAffectedNodes4");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = *(Suitors.at(x)->partners.end() - Suitors.at(x)->num_visits);
//         Suitors.at(x)->num_visits--;
//         // const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//         INFO("Inserting suitor ", x, " into ", y.id);
//         auto weight = G->weight(y.id, x);
//         auto minWeight = Suitors.at(y.id)->min.weight;
//         if (minWeight < weight || (minWeight == weight && Suitors.at(y.id)->min.id > x)) {
//             Proposed.at(x)->insert({y.id, y.weight});
//             DynBNode newS = Suitors.at(y.id)->insert({x, y.weight});

//             if (newS.id != none) {
//                 Suitors.at(newS.id)->remove(y.id);
//                 Proposed.at(y.id)->remove(newS.id);
//             }
//         }
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

// void DynamicBSuitorMatcher::updateAffectedNodes5() {
//     affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
//     INFO("Starting updateAffectedNodes5");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//         INFO("Inserting suitor ", x, " into ", y.id);
//         auto weight = G->weight(y.id, x);
//         auto minWeight = Suitors.at(y.id)->min.weight;
//         if (minWeight < weight || (minWeight == weight && Suitors.at(y.id)->min.id > x)) {
//             DynBNode newS = Suitors.at(y.id)->insert({x, y.weight});
//             if (newS.id != none) {
//                 Suitors.at(newS.id)->remove(y.id);
//             }
//         }
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

// void DynamicBSuitorMatcher::updateAffectedNodes6() {
//     affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
//     INFO("Starting updateAffectedNodes6");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//         INFO("Inserting suitor ", x, " into ", y.id);
//         auto weight = G->weight(y.id, x);
//         auto minWeight = Suitors.at(y.id)->min.weight;
//         if (minWeight < weight || (minWeight == weight && Suitors.at(y.id)->min.id > x)) {
//             DynBNode newS = Suitors.at(y.id)->insert({x, y.weight});
//             if (newS.id != none) {
//                 INFO("Updating loose end ", newS.id);
//                 Suitors.at(newS.id)->remove(y.id);

//                 DynBNode looseEnd = DynBNode{none, Suitors.at(newS.id)->min.weight};

//                 G->forNeighborsOf(newS.id, [&](node x, edgeweight edgeWeight) {
//                     // ignore edges that that still need to be processed
//                     if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(newS.id, x))
//                             != edgeBatch.end()
//                         || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, newS.id))
//                                != edgeBatch.end())
//                         return;
//                     if (Suitors.at(newS.id)->hasPartner(x))
//                         return;

//                     const auto z = Suitors.at(x)->min;

//                     if (z.id == none && looseEnd.weight < edgeWeight) {
//                         looseEnd.id = x;
//                         looseEnd.weight = edgeWeight;
//                     }
//                 });

//                 if (looseEnd.id != none) {
//                     INFO("Found a loose end ", looseEnd.id, " with weight ", looseEnd.weight);
//                     Suitors.at(newS.id)->insert({looseEnd.id, looseEnd.weight});
//                     Suitors.at(looseEnd.id)->insert({newS.id, looseEnd.weight});
//                 }
//             }
//         }
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

// void DynamicBSuitorMatcher::updateAffectedNodes7(Operation op) {
//     affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
//     INFO("Starting updateAffectedNodes7");
//     while (affectedNodes.size() > 1) {
//         const node x = affectedNodes.back();
//         const auto y = Suitors.at(x)->partners.back();
//         affectedNodes.pop_back();
//         INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
//         INFO("Inserting suitor ", x, " into ", y.id);
//         auto weight = G->weight(y.id, x);
//         auto minWeight = Suitors.at(y.id)->min.weight;
//         if (minWeight < weight || (minWeight == weight && Suitors.at(y.id)->min.id > x)) {
//             DynBNode newS = Suitors.at(y.id)->insert({x, y.weight});
//             if (newS.id != none && op == Operation::Insert) {
//                 INFO("Updating loose end ", newS.id);
//                 Suitors.at(newS.id)->remove(y.id);

//                 DynBNode looseEnd = DynBNode{none, Suitors.at(newS.id)->min.weight};

//                 G->forNeighborsOf(newS.id, [&](node x, edgeweight edgeWeight) {
//                     // ignore edges that that still need to be processed
//                     if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(newS.id, x))
//                             != edgeBatch.end()
//                         || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, newS.id))
//                                != edgeBatch.end())
//                         return;
//                     if (Suitors.at(newS.id)->hasPartner(x))
//                         return;

//                     const auto z = Suitors.at(x)->min;

//                     if (z.id == none && looseEnd.weight < edgeWeight) {
//                         looseEnd.id = x;
//                         looseEnd.weight = edgeWeight;
//                     }
//                 });

//                 if (looseEnd.id != none) {
//                     INFO("Found a loose end ", looseEnd.id, " with weight ", looseEnd.weight);
//                     Suitors.at(newS.id)->insert({looseEnd.id, looseEnd.weight});
//                     Suitors.at(looseEnd.id)->insert({newS.id, looseEnd.weight});
//                 }
//             }
//         }
//         affected[y.id] = false;
//         affected[x] = false;
//     }
// }

void DynamicBSuitorMatcher::updateAffectedNodes(uint8_t batchId) {
    affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
    INFO("Starting updateAffectedNodes8");
    while (affectedNodes.size() > 0) {
        const node x = affectedNodes.back();
        const auto y = Suitors.at(x)->popLastUpdate();
        affectedNodes.pop_back();
        INFO("Trying to add ", x, " to ", y.id);

        if (y.id == none || !Suitors.at(x)->hasPartner(y.id)) {
            INFO("Either none or ", y.id, " not in ", x, " anymore.");
            continue;
        }
        auto weight = G->weight(y.id, x);
        auto minWeight = Suitors.at(y.id)->min.weight;

        if (minWeight > weight || (minWeight == weight && Suitors.at(y.id)->min.id < x)) {
            INFO("Not adding ", x, " to ", y.id, ". Searching for new partner for ", x);
            Suitors.at(x)->remove(y.id);

            DynBNode looseEnd = DynBNode{none, Suitors.at(x)->min.weight};

            bool done = false;
            auto current = x;
            edgeweight prev = std::numeric_limits<edgeweight>::max();

            do {

                done = true;

                // Check for new candidate
                G->forNeighborsOf(current, [&](node v, edgeweight edgeWeight) {
                    auto edgeHash = current < v ? std::make_pair(current, v) : std::make_pair(v, current);
                    // ignore edges that that still need to be processed
                    if (batchTracker.contains(edgeHash)) {
                        if (batchId < batchTracker[edgeHash])
                            return;
                    }

                    // if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, v))
                    //         != edgeBatch.end()
                    //     || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(v, current))
                    //            != edgeBatch.end())
                    //     return;
                    if (Suitors.at(current)->hasPartner(v))
                        return;

                    const auto z = Suitors.at(v)->min;

                    INFO("Node: ", current, "(weight:", looseEnd.weight, ") Evaluating ", v,
                         " (edge_weight: ", edgeWeight, " min: ", z.id, ", weight: ", z.weight, ")");
                    if ((edgeWeight > looseEnd.weight
                         || (looseEnd.weight == edgeWeight && v < looseEnd.id))
                        && (edgeWeight > z.weight || (edgeWeight == z.weight && current < z.id))
                        && (edgeWeight <= prev)) {
                        looseEnd.id = v;
                        looseEnd.weight = edgeWeight;
                        return;
                    }
                });

                if (looseEnd.id == none || looseEnd.weight < Suitors.at(looseEnd.id)->min.weight
                    || (looseEnd.weight == Suitors.at(looseEnd.id)->min.weight
                        && current > Suitors.at(looseEnd.id)->min.id)) {
                    break;
                }

                INFO("Done searching neighbors: \n cur: ", current, ", partner: ", looseEnd.id,
                     ", weight: ", looseEnd.weight);
                DynBNode problem =
                    Suitors.at(current)->insert({looseEnd.id, looseEnd.weight}, false);
                if (problem.id != none) {
                    INFO("Pushed ", looseEnd.id, " to ", current, " but it was full with min ",
                         problem.id);
                }
                INFO("Suitors at y ", looseEnd.id, " before update: ", *Suitors.at(looseEnd.id));
                DynBNode y = Suitors.at(looseEnd.id)->insert({current, looseEnd.weight}, false);
                INFO("Suitors at y ", looseEnd.id, " after update: ", *Suitors.at(looseEnd.id));

                // line 16-20
                if (y.id != none) {
                    INFO("Removing suitor ", looseEnd.id, " from ", y.id, " (new loose end).");
                    Suitors.at(y.id)->remove(looseEnd.id);
                    current = y.id;
                    done = false;
                }

                // line 21-23
                prev = looseEnd.weight;
                looseEnd.id = Suitors.at(current)->min.id;
                looseEnd.weight = Suitors.at(current)->min.weight;

            } while (!done);

            // affected[y.id] = false;
            // affected[x] = false;
            continue;
        }

        INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
        INFO("Inserting suitor ", x, " into ", y.id);
        DynBNode newS = Suitors.at(y.id)->insert({x, weight}, false);
        if (newS.id != none) {
            INFO("Updating loose end ", newS.id, " of ", y.id);
            Suitors.at(newS.id)->remove(y.id);

            DynBNode looseEnd = DynBNode{none, Suitors.at(newS.id)->min.weight};

            bool done = false;
            auto current = newS.id;
            edgeweight prev = std::numeric_limits<edgeweight>::max();

            // Create a path for loose ends
            // No need here to cover affected nodes, since we don't have two edge ends
            // In each step the s-variant is maintained
            do {

                done = true;

                // Check for new candidate
                G->forNeighborsOf(current, [&](node v, edgeweight edgeWeight) {
                    auto edgeHash = current < v ? std::make_pair(current, v) : std::make_pair(v, current);
                    // ignore edges that that still need to be processed
                    if (batchTracker.contains(edgeHash)) {
                        if (batchId < batchTracker[edgeHash])
                            return;
                    }                    // if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, v))
                    //         != edgeBatch.end()
                    //     || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(v, current))
                    //            != edgeBatch.end())
                    //     return;
                    if (Suitors.at(current)->hasPartner(v))
                        return;

                    const auto z = Suitors.at(v)->min;

                    INFO("Node: ", current, "(weight:", looseEnd.weight, ") Evaluating ", v,
                         " (edge_weight: ", edgeWeight, " min: ", z.id, ", weight: ", z.weight, ")");
                    if ((edgeWeight > looseEnd.weight
                         || (looseEnd.weight == edgeWeight && v < looseEnd.id))
                        && (edgeWeight > z.weight || (edgeWeight == z.weight && current < z.id))
                        && (edgeWeight <= prev)) {
                        looseEnd.id = v;
                        looseEnd.weight = edgeWeight;
                        return;
                    }

                });

                if (looseEnd.id == none || looseEnd.weight < Suitors.at(looseEnd.id)->min.weight
                    || (looseEnd.weight == Suitors.at(looseEnd.id)->min.weight
                        && current > Suitors.at(looseEnd.id)->min.id)) {
                    break;
                }

                INFO("Done searching neighbors: \n cur: ", current, ", partner: ", looseEnd.id,
                     ", weight: ", looseEnd.weight);
                DynBNode problem =
                    Suitors.at(current)->insert({looseEnd.id, looseEnd.weight}, false);
                if (problem.id != none) {
                    INFO("Pushed ", looseEnd.id, " to ", current, " but it was full with min ",
                         problem.id);
                }
                INFO("Suitors at y ", looseEnd.id, " before update: ", *Suitors.at(looseEnd.id));
                DynBNode y = Suitors.at(looseEnd.id)->insert({current, looseEnd.weight}, false);
                INFO("Suitors at y ", looseEnd.id, " after update: ", *Suitors.at(looseEnd.id));

                // line 16-20
                if (y.id != none) {
                    INFO("Removing suitor ", looseEnd.id, " from ", y.id, " (new loose end).");
                    Suitors.at(y.id)->remove(looseEnd.id);
                    current = y.id;
                    done = false;
                }

                // line 21-23
                prev = looseEnd.weight;
                looseEnd.id = Suitors.at(current)->min.id;
                looseEnd.weight = Suitors.at(current)->min.weight;

            } while (!done);
        }
        // }
        // affected[y.id] = false;
        // affected[x] = false;
    }
}

count DynamicBSuitorMatcher::getNumberOfAffected() {
    return affectedNodesPerRun;
}

} // namespace NetworKit
