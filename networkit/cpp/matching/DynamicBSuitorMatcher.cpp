#include <cassert>
#include <chrono>

#include <networkit/matching/DynamicBSuitorMatcher.hpp>

namespace NetworKit {

void DynamicBSuitorMatcher::processEdgeInsertion(const WeightedEdge &edge) {

    node u = edge.u;
    node v = edge.v;
    edgeweight w = G->weight(u,v);
    
    DynBNode startU = Suitors.at(u)->insert({v,w});
    DynBNode startV = Suitors.at(v)->insert({u,w});
    affectedNodesPerRun += 2;
    INFO("StartU: ", startU.id);
    INFO("StartV: ", startV.id);

    if (startU.id != none) {
        Suitors.at(startU.id)->remove(u);
    }

    if (startV.id != none) {
        Suitors.at(startV.id)->remove(v);
    }

    if(startU.id != none) {
        trackUpdatePath(0, startU.id);
    }

    if(startV.id != none) {
        trackUpdatePath(0, startV.id);
    }
}

void DynamicBSuitorMatcher::trackUpdatePath(size_t batchId, node start, bool recursiveCall) {
    INFO("Start node: ", start);

    bool done = false;

    node current = start;
    node partner = Suitors.at(current)->min.id;
    auto heaviest = Suitors.at(current)->min.weight;
    INFO("Starting findAffectedNodes8 for id ", batchId, " with: \n cur: ", current, ", partner: ", partner,
         ", weight: ", heaviest);

    edgeweight prev = std::numeric_limits<edgeweight>::max();

    std::vector<DynBNode> looseEnds;


    do {
        done = true;

        G->forNeighborsOf(current, [&](node x, edgeweight weight) {
            if (Suitors.at(current)->hasPartner(x))
                return;

            const auto z = Suitors.at(x)->min;
            INFO("Node: ", current, " (weight:", heaviest, ") Evaluating ", x, " (edge_weight: ", weight, " min: ", z.id,
                 ", min_weight: ", z.weight, ")");

            if ((weight > heaviest || (weight == heaviest && x < partner))
                && (weight > z.weight || (weight == z.weight && current < z.id))
                && (weight <= prev)) {
                partner = x;
                heaviest = weight;
                return;
            }
        });

        if (partner == none || heaviest < Suitors.at(partner)->min.weight
            || (heaviest == Suitors.at(partner)->min.weight
                && current > Suitors.at(partner)->min.id)) {
            break;
        }

        INFO("Suitors at current ", current, ": ", *Suitors.at(current));
        INFO("Suitors at partner ", partner, ": ", *Suitors.at(partner));
        INFO("Min at partner ", partner, ": ", (Suitors.at(partner)->min).id,
             " (weight: ", (Suitors.at(partner)->min).weight, ")");

        INFO("Inserting suitor ", current, " into ", partner, " and vice versa.");
        DynBNode prevCurrent = Suitors.at(current)->insert({partner, heaviest});
        DynBNode prevPartner = Suitors.at(partner)->insert({current, heaviest});
        affectedNodesPerRun++;

        if(prevCurrent.id != none) {
            INFO("Current was saturated. Removing current from prevCurrent ", prevCurrent.id);
            Suitors.at(prevCurrent.id)->remove(current);
            looseEnds.emplace_back(DynBNode{prevCurrent.id,Suitors.at(prevCurrent.id)->min.weight});
        }


        if (prevPartner.id != none) {
            INFO("Removing suitor ", prevPartner.id, " from ", partner);
            INFO("Suitors at y ", prevPartner.id, ": ", *Suitors.at(prevPartner.id));

            INFO("Removing suitor ", partner, " from ", prevPartner.id);
            Suitors.at(prevPartner.id)->remove(partner);
            current = prevPartner.id;
            done = false;
        }

        // line 21-23
        prev = heaviest;
        partner = Suitors.at(current)->min.id;
        heaviest = Suitors.at(current)->min.weight;
    } while (!done);

    for (auto &looseEnd : looseEnds) {
        
        trackUpdatePath(batchId, looseEnd.id, true);
    }
}

void DynamicBSuitorMatcher::processEdgeRemoval(const Edge &edge) {

    node u = edge.u;
    node v = edge.v;
    
    Suitors.at(u)->remove(v);
    Suitors.at(v)->remove(u);

    trackUpdatePath(0, u);    
    trackUpdatePath(0, v);
}


void DynamicBSuitorMatcher::addEdges(std::vector<WeightedEdge> &edges, bool sort) {

    affectedNodesPerRun = 0;

    if (sort) {
        std::sort(edges.begin(), edges.end(),
                [](const WeightedEdge &a, const WeightedEdge &b) { return a.weight > b.weight; });
    }

    for (const auto &edge : edges) {
        if ((Suitors.at(edge.u)->hasPartner(edge.v) && Suitors.at(edge.v)->hasPartner(edge.u))
            || !isBetterMatch(edge.u, edge.v, edge.weight)
            || !isBetterMatch(edge.v, edge.u, edge.weight)) {
            INFO("Edge ", edge.u, ",", edge.v,
                 " ignored, since min(u): ", Suitors.at(edge.u)->min.weight,
                 " min(v): ", Suitors.at(edge.v)->min.weight);
            continue;
        }
        INFO("Edge ", edge.u, ",", edge.v, " better choice. Min at ", edge.u, ": ", Suitors.at(edge.u)->min.weight, " Min at ", edge.v, ": ", Suitors.at(edge.v)->min.weight);

        processEdgeInsertion(edge);
    }
}

void DynamicBSuitorMatcher::removeEdges(std::vector<Edge> &edges) {

    affectedNodesPerRun = 0;

    for (const auto &edge : edges) {
        assert(!G->hasEdge(edge.u, edge.v));
        if (Suitors.at(edge.u)->hasPartner(edge.v)) {

            processEdgeRemoval(edge);
        }
    }
}

count DynamicBSuitorMatcher::getNumberOfAffected() {
    return affectedNodesPerRun;
}

} // namespace NetworKit
