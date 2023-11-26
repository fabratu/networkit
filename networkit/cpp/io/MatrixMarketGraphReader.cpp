#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/MatrixMarketGraphReader.hpp>
#include <networkit/io/MatrixMarketParser.hpp>

namespace NetworKit {

Graph MatrixMarketGraphReader::read(const std::string &path) {
    MatrixMarketParser parser(path);

    auto header = parser.getHeader();
    auto size = parser.getMatrixSize();
    auto weighted = true;
    auto directed = false;

    if (header.field == MatrixMarketParser::Field::Pattern)
        weighted = false;
    if (header.symmetry == MatrixMarketParser::Symmetry::General)
        directed = true;

    Graph G(std::max(size.rows, size.columns), weighted, directed);
    std::string graphName =
        Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

    std::optional<MatrixMarketParser::Edge> current_edge = parser.getNext(weighted);
    count edgeCount = 0;
    count selfLoopCount = 0;

    while (current_edge.has_value()) {
        const auto [from, to, weight] = current_edge.value();
        edgeCount++;
        if (from == to)
            selfLoopCount++;

        if (weighted && directed) {
            G.addEdge(from, to, *weight);
        } else if (weighted && !directed) {
            G.addPartialEdge(unsafe, from, to, *weight);
            G.addPartialEdge(unsafe, to, from, *weight);
        } else if (!weighted && directed) {
            G.addEdge(from, to);
        } else {
            G.addPartialEdge(unsafe, from, to);
            G.addPartialEdge(unsafe, to, from);
        }

        current_edge = parser.getNext(weighted);
    }

    G.setEdgeCount(unsafe, edgeCount);
    G.setNumberOfSelfLoops(unsafe, selfLoopCount);

    if (header.format == MatrixMarketParser::Format::Coordinate && G.numberOfEdges() != size.nonzeros) {
        ERROR("Number of added edges doesn't match the specifed number of edges in file ", path);
    }

    return G;
}

} /* namespace NetworKit */
