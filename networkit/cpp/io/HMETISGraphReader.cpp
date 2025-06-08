/*
 * METISGraphReader.cpp
 *
 *  Created on: 05.06.2025
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/auxiliary/NumberParsing.hpp>
#include <networkit/io/HMETISGraphReader.hpp>

namespace NetworKit {

HMETISGraphReader::HMETISGraphReader(node firstNode) : offset(firstNode) {}

std::vector<node> HMETISGraphReader::parseLine(std::string_view line, count ignoreFirst) {
    auto it = line.begin();
    auto end = line.end();
    std::vector<node> adjacencies;
    node v;
    index i = 0;

    while (i < ignoreFirst) {
        // parse first values but ignore them.
        double dummy;
        std::tie(dummy, it) = Aux::Parsing::strTo<double>(it, end);
        ++i;
    }

    while (it != end) {
        std::tie(v, it) = Aux::Parsing::strTo<node>(it, end);
        // According to https://github.com/Coloquinte/minipart/
        // at least some instances (e.g. ibm-curcuits) are 1-indexed.
        // However, the hmetis-format description does not support this
        // assumption. Therefore 0-indexed is the default case.
        adjacencies.push_back(v - offset);
    }

    return adjacencies;
}

Hypergraph HMETISGraphReader::read(std::string_view path) {
    hGraphFile.open(path.data());
    if (!hGraphFile.is_open()) {
        throw std::runtime_error("Could not open file: " + std::string(path));
    }

    // Read the header
    std::string line;
    if (!std::getline(hGraphFile, line)) {
        throw std::runtime_error("Failed to read header from file: " + std::string(path));
    }

    auto parts = Aux::StringTools::split(line, ' ');
    if (parts.size() < 2) {
        throw std::runtime_error("Invalid header format in file: " + std::string(path));
    }

    count n = std::stoul(parts[0]);
    count m = std::stoul(parts[1]);

    Hypergraph hGraph(n, 0);

    // Read the hyperedges
    while (std::getline(hGraphFile, line)) {
        hGraph.addEdge(parseLine(line, 0));

        // auto edgeNodes = Aux::StringTools::split(line, ' ');
        // std::vector<node> edge;
        // for (const auto &nodeStr : edgeNodes) {
        //     edge.push_back(std::stoul(nodeStr));
        // }
        // hGraph.addEdge(edge);
    }

    hGraphFile.close();
    return hGraph;
}

} // namespace NetworKit
