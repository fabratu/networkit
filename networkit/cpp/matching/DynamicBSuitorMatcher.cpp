#include <cassert>

#include <networkit/matching/DynamicBSuitorMatcher.hpp>

namespace NetworKit {

void DynamicBSuitorMatcher::processEdgeInsertion(const WeightedEdge &edge) {
    affected[edge.u] = true;
    affectedNodes.clear();
    findAffectedNodes(edge.u, edge.v, Operation::Insert);
    updateAffectedNodes();
    affected[edge.u] = false;

    affected[edge.v] = true;
    affectedNodes.clear();
    findAffectedNodes(edge.v, edge.u, Operation::Insert);
    updateAffectedNodes();

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

void DynamicBSuitorMatcher::findAffectedNodes(node u, node v, Operation op) {
    bool done = false;

    node current = u;
    node partner = (op == Operation::Insert) ? v : none;
    auto heaviest = G->weight(current, v);
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
            if (Suitors.at(current)->hasPartner(x))
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

        INFO("Suitors at current ", current, ": ", *Suitors.at(current));
        INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
        INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
             " (weight: ", (Suitors.at(partner)->min).weight, ")");

        const auto ps = Suitors.at(partner)->min;
        if (heaviest > ps.weight || (heaviest == ps.weight && current < ps.id)) {
            const auto y = Suitors.at(partner)->popMinIfFull();

            if (y.id != none) {
                INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
                INFO("Removing suitor ", y.id, " from ", partner);
            }

            INFO("Inserting suitor ", current, " into ", partner);
            Node newS = Suitors.at(partner)->insert({current, heaviest});
            if (newS.id != none) {
                Suitors.at(newS.id)->remove(partner);
            }
            affected[partner] = true;
            affectedNodes.emplace_back(partner);

            if (y.id != none) {
                INFO("Removing suitor ", partner, " from ", y.id);
                Suitors.at(y.id)->remove(partner);
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
    while (affectedNodes.size() > 1) {
        const node x = affectedNodes.back();
        const auto y = Suitors.at(x)->partners.back();
        affectedNodes.pop_back();
        INFO("Suitors at y ", y.id, ": ", *Suitors.at(y.id));
        INFO("Inserting suitor ", x, " into ", y.id);
        Node newS = Suitors.at(y.id)->insert({x, y.weight});
        if (newS.id != none) {
            Suitors.at(newS.id)->remove(y.id);
        }
        affected[y.id] = false;
        affected[x] = false;
    }
}

count DynamicBSuitorMatcher::getNumberOfAffected() {
    return affectedNodesPerRun.size();
}

} // namespace NetworKit
