/*
 * HMETISHypergraphReader.cpp
 */

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/HMETISHypergraphReader.hpp>

namespace NetworKit {
namespace {

bool isComment(const std::string &line) {
    const auto first = line.find_first_not_of(" \t\r");
    return first != std::string::npos && line[first] == '%';
}

bool isBlank(const std::string &line) {
    return line.find_first_not_of(" \t\r") == std::string::npos;
}

std::runtime_error parseError(std::string_view path, count line, const std::string &message) {
    return std::runtime_error{"Invalid hMETIS file '" + std::string{path} + "' at line "
                              + std::to_string(line) + ": " + message};
}

} // namespace

Hypergraph HMETISHypergraphReader::read(std::string_view path) const {
    std::ifstream file{std::string{path}};
    Aux::enforceOpened(file);

    std::string line;
    count lineNumber = 0;
    while (std::getline(file, line)) {
        ++lineNumber;
        if (!isComment(line) && !isBlank(line))
            break;
    }
    if (isComment(line) || isBlank(line))
        throw parseError(path, lineNumber, "missing header");

    count numberOfEdges;
    count numberOfNodes;
    index format = 0;
    std::string trailing;
    std::istringstream header{line};
    if (!(header >> numberOfEdges >> numberOfNodes))
        throw parseError(path, lineNumber, "expected '<hyperedges> <nodes> [format]'");
    if (header >> format) {
        if (header >> trailing)
            throw parseError(path, lineNumber, "too many header fields");
    } else if (!header.eof()) {
        throw parseError(path, lineNumber, "invalid format flag");
    }

    if (format != 0 && format != 1 && format != 20 && format != 21)
        throw parseError(path, lineNumber, "unsupported format flag (expected 0, 1, 20, or 21)");

    const bool edgeWeighted = format % 10 == 1;
    const bool incidenceWeighted = format / 10 == 2;
    Hypergraph hypergraph{numberOfNodes, 0, edgeWeighted, incidenceWeighted};

    count edgesRead = 0;
    while (edgesRead < numberOfEdges && std::getline(file, line)) {
        ++lineNumber;
        if (isComment(line))
            continue;

        std::istringstream edgeLine{line};
        edgeweight edgeWeight = defaultEdgeWeight;
        if (edgeWeighted && !(edgeLine >> edgeWeight))
            throw parseError(path, lineNumber, "missing hyperedge weight");

        std::vector<node> members;
        std::vector<edgeweight> weights;
        std::unordered_set<node> uniqueMembers;
        count fileNode;
        while (edgeLine >> fileNode) {
            if (fileNode == 0 || fileNode > numberOfNodes)
                throw parseError(path, lineNumber, "node id is outside the declared range");

            const node u = fileNode - 1;
            if (!uniqueMembers.insert(u).second)
                throw parseError(path, lineNumber, "duplicate node id in a hyperedge");

            members.push_back(u);
            if (incidenceWeighted) {
                edgeweight weight;
                if (!(edgeLine >> weight))
                    throw parseError(path, lineNumber, "missing incidence weight");
                weights.push_back(weight);
            }
        }
        if (!edgeLine.eof())
            throw parseError(path, lineNumber, "invalid hyperedge field");

        const edgeid eid = hypergraph.addEdge(members);
        if (edgeWeighted)
            hypergraph.setEdgeWeight(eid, edgeWeight);
        if (incidenceWeighted) {
            for (index i = 0; i < members.size(); ++i)
                hypergraph.setIncidenceWeight(members[i], eid, weights[i]);
        }
        ++edgesRead;
    }

    if (edgesRead != numberOfEdges)
        throw parseError(path, lineNumber, "fewer hyperedges than declared in the header");

    while (std::getline(file, line)) {
        ++lineNumber;
        if (!isComment(line) && !isBlank(line))
            throw parseError(path, lineNumber, "more hyperedges than declared in the header");
    }

    return hypergraph;
}

} // namespace NetworKit
