/*
 * HMETISGraphReader.hpp
 *
 *  Created on: 05.06.2025
 *      Author: Fabian Brandt-Tumescheit
 *
 *  Code partially reused from https://github.com/Coloquinte/minipart
 */

#ifndef NETWORKIT_IO_HMETIS_GRAPH_READER_HPP_
#define NETWORKIT_IO_HMETIS_GRAPH_READER_HPP_

#include <networkit/io/HypergraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the hMETIS file format documented in [1]
 *
 * [1] https://course.ece.cmu.edu/~ee760/760docs/hMetisManual.pdf
 */
class HMETISGraphReader final : public HypergraphReader {
public:
    HMETISGraphReader() = default;

    /**
     * @param[in]  firstNode  index of the first node in the file
     */
    HMETISGraphReader(node firstNode);

    /**
     * Takes a file path as parameter and returns a graph file.
     *
     * @param[in]  path  file path
     *
     * @param[out]  the graph read from file
     */
    Hypergraph read(std::string_view path) override;

private:
    std::vector<node> parseLine(std::string_view line, count ignoreFirst = 0);
    std::ifstream hGraphFile;
    node offset = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_HMETIS_GRAPH_READER_HPP_
