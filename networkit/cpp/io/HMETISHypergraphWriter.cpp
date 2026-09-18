/*
 * HMETISHypergraphWriter.cpp
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/HMETISHypergraphWriter.hpp>

namespace NetworKit {

void HMETISHypergraphWriter::write(const Hypergraph &hypergraph, std::string_view path) const {
    std::ofstream file{std::string{path}};
    Aux::enforceOpened(file);
    file << std::setprecision(std::numeric_limits<edgeweight>::max_digits10);

    std::unordered_map<node, node> nodeIds;
    node nextNode = 1;
    hypergraph.forNodes([&](node u) { nodeIds.emplace(u, nextNode++); });

    const index format = static_cast<index>(hypergraph.isWeighted())
                         + 20 * static_cast<index>(hypergraph.isIncidenceWeighted());
    file << hypergraph.numberOfEdges() << ' ' << hypergraph.numberOfNodes() << ' ' << format
         << '\n';

    hypergraph.forEdges([&](edgeid eid) {
        if (hypergraph.isWeighted())
            file << hypergraph.getEdgeWeight(eid);

        std::vector<node> members{hypergraph.edgeMembers(eid).begin(),
                                  hypergraph.edgeMembers(eid).end()};
        std::sort(members.begin(), members.end(),
                  [&](node u, node v) { return nodeIds.at(u) < nodeIds.at(v); });

        for (node u : members) {
            if (hypergraph.isWeighted() || u != members.front())
                file << ' ';
            file << nodeIds.at(u);
            if (hypergraph.isIncidenceWeighted())
                file << ' ' << hypergraph.getIncidenceWeight(u, eid);
        }
        file << '\n';
    });
}

} // namespace NetworKit
