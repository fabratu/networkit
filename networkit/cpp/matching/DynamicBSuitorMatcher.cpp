#include <cassert>

#include <networkit/matching/DynamicBSuitorMatcher.hpp>

namespace NetworKit {

void DynamicBSuitorMatcher::processEdgeInsertion(const WeightedEdge &edge) {
    affected[edge.u] = true;
    affectedNodes.clear();
    findAffectedNodes(edge.u, edge.v, Operation::Insert, true);
    // updateAffectedNodes();
    affected[edge.u] = false;

    affected[edge.v] = true;
    affectedNodes.clear();
    findAffectedNodes(edge.v, edge.u, Operation::Insert);
    // updateAffectedNodes();

    affected[edge.u] = affected[edge.v] = false;
}

void DynamicBSuitorMatcher::processEdgeRemoval(const Edge &edge) {
    Suitors.at(edge.u)->remove(edge.v);
    Suitors.at(edge.v)->remove(edge.u);

    affected[edge.u] = affected[edge.v] = true;
    affectedNodes.clear();
    affectedNodes.emplace_back(edge.u);
    findAffectedNodes(edge.u, edge.v, Operation::Remove);
    updateAffectedNodes();
    affected[edge.u] = false;

    affected[edge.v] = true;
    affectedNodes.clear();
    affectedNodes.emplace_back(edge.v);
    findAffectedNodes(edge.v, edge.u, Operation::Remove);
    updateAffectedNodes();
    affected[edge.v] = false;
}

void DynamicBSuitorMatcher::addEdges(std::vector<WeightedEdge> &edges) {
    auto convertToEdge = [](const WeightedEdge &weightedEdge) {
        return Edge(weightedEdge.u, weightedEdge.v);
    };
    std::transform(edges.begin(), edges.end(), std::back_inserter(edgeBatch), convertToEdge);
    std::sort(edges.begin(), edges.end(),
              [](const WeightedEdge &a, const WeightedEdge &b) { return a.weight > b.weight; });

    for (const auto &edge : edges) {
        if ((Suitors.at(edge.u)->hasPartner(edge.v) && Suitors.at(edge.v)->hasPartner(edge.u))
            || !isBetterMatch(edge.u, edge.v, edge.weight)
            || !isBetterMatch(edge.v, edge.u, edge.weight))
            continue;

        affectedNodes.clear();

        processEdgeInsertion(edge);
        affected[edge.u] = affected[edge.v] = false;

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
    edgeBatch = edges;

    for (const auto &edge : edges) {
        assert(!G->hasEdge(edge.u, edge.v));
        if (Suitors.at(edge.u)->hasPartner(edge.v)) {

            affectedNodes.clear();
            assert(numberOfAffectedEquals(0));

            processEdgeRemoval(edge);
            affected[edge.u] = affected[edge.v] = false;

            G->forNodes([&](node u) {
                for (auto s : Suitors.at(u)->partners) {
                    assert(Suitors.at(s.id)->hasPartner(u));
                }
            });
        }
    }
}

void DynamicBSuitorMatcher::findAffectedNodes(node u, node v, Operation op, bool skipFirst) {
    bool done = false;
    bool localSkip = skipFirst;
    node current = u;
    node partner = (op == Operation::Insert) ? v : none;
    auto heaviest = G->weight(current, v);
    node toRemoveIfRevisited = none;
    INFO("Starting findAffectedNodes with: \n cur: ", current, ", partner: ", partner,
         ", weight: ", heaviest);

    edgeweight prev = std::numeric_limits<edgeweight>::max();

    do {
        G->forNeighborsOf(current, [&](node x, edgeweight weight) {
            // ignore edges that that still need to be processed
            if (std::find(edgeBatch.begin(), edgeBatch.end(), Edge(current, x)) != edgeBatch.end()
                || std::find(edgeBatch.begin(), edgeBatch.end(), Edge(x, current))
                       != edgeBatch.end())
                return;
            if (Proposed.at(current)->hasPartner(x))
                return;

            const auto z = Suitors.at(x)->min;
            if (!affected[x] || (op == Operation::Insert && (x == u || x == v))) {
                if ((weight > heaviest || (weight == heaviest && x < partner))
                    && (weight > z.weight || (weight == z.weight && current < z.id))
                    && (weight <= prev)) {
                    partner = x;
                    heaviest = weight;
                    return;
                }
            }
        });
        done = true;
        INFO("Done searching neighbors: \n cur: ", current, ", partner: ", partner,
             ", weight: ", heaviest);

        if (partner == none) {
            affected[current] = false;
            done = true;
            continue;
        }
        INFO("Proposed at partner ", partner, ": ", *Proposed.at(partner));
        INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
        INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
             " (weight: ", (Suitors.at(partner)->min).weight, ")");
        const auto ps = Suitors.at(partner)->min;
        if (heaviest > ps.weight || (heaviest == ps.weight && current < ps.id)) {
            const auto y = Suitors.at(partner)->popMinIfFull();
            INFO("Removing suitor ", y.id, " from ", partner);
            if (y.id != none) {
                // INFO("Removing proposed ", partner, " from ", y.id);
                // Proposed.at(y.id)->remove(partner);
                INFO("Removing suitor ", partner, " from ", y.id);
                Suitors.at(y.id)->remove(partner);
                INFO("Removing proposed ", y.id, " from ", partner);
                Proposed.at(partner)->remove(y.id);
            }

            // const auto cs = Suitors.at(current)->getMinAnyway();
            // if (cs.id != none && Proposed.at(current)->isFull()) {
            //     INFO("[CURRENT] Removing proposed ", cs.id, " from ", current);
            //     Proposed.at(current)->remove(cs.id);
            // }

            INFO("Inserting suitor ", current, " into ", partner, " with weight ", heaviest);
            Suitors.at(partner)->insert({current, heaviest});
            INFO("Inserting proposed ", partner, " into ", current);
            if (current == u)
                Proposed.at(current)->insert({partner, heaviest}, true);
            else
                Proposed.at(current)->insert({partner, heaviest});
            affected[partner] = true;
            affectedNodes.emplace_back(partner);

            INFO("Inserting proposed ", current, " into ", partner);
            Proposed.at(partner)->insert({current, heaviest});
            if (!localSkip) {
                INFO("Inserting suitor ", partner, " into ", current, " with weight ", heaviest);
                Suitors.at(current)->insert({partner, heaviest});
                localSkip = false;
            }

            // if (y.id == none) {
            //     INFO("Inserting proposed ", current, " into ", partner);
            //     Proposed.at(partner)->insert({current, heaviest});
            //     INFO("Inserting suitor ", partner, " into ", current, " with weight ", heaviest);
            //     Suitors.at(current)->insert({partner, heaviest});
            // }

            if (y.id != none) {
                affected[y.id] = true;
                current = y.id;
                done = false;
            }
        } else {
            affected[current] = false;
        }
        prev = heaviest;
        partner = Suitors.at(current)->min.id;
        heaviest = Suitors.at(current)->min.weight;
        INFO("End of current loop: \n cur: ", current, ", partner: ", partner,
             ", weight: ", heaviest);
    } while (!done);
}

void DynamicBSuitorMatcher::updateAffectedNodes() {
    affectedNodesPerRun.insert(affectedNodes.begin(), affectedNodes.end());
    INFO("Starting updateAffectedNodes");
    for (node x : affectedNodes) {
        const auto y = Suitors.at(x)->partners.back();
        INFO("Inserting suitor ", x, " into ", y.id);
        Suitors.at(y.id)->insert({x, y.weight});
        INFO("Inserting proposed ", y.id, " into ", x);
        Proposed.at(x)->insert({y.id, y.weight});
        affected[y.id] = false;
        affected[x] = false;
    }
    affectedNodes.clear();
    // while (affectedNodes.size() > 1) {
    //     const node x = affectedNodes.front();
    //     const auto y = Suitors.at(x)->partners.back();
    //     affectedNodes.pop_back();
    //     INFO("Inserting suitor ", x, " into ", y.id);
    //     Suitors.at(y.id)->insert({x, y.weight});
    //     INFO("Inserting proposed ", y.id, " into ", x);
    //     Proposed.at(x)->insert({y.id, y.weight});
    //     affected[y.id] = false;
    //     affected[x] = false;
    // }
}

count DynamicBSuitorMatcher::getNumberOfAffected() {
    return affectedNodesPerRun.size();
}

} // namespace NetworKit
