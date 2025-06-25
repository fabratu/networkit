/*
 * GraphReader.hpp
 *
 *  Created on: 05.06.2025
 *      Author: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_IO_HYPERGRAPH_READER_HPP_
#define NETWORKIT_IO_HYPERGRAPH_READER_HPP_

#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {
/**
 * @ingroup io
 * Abstract base class for graph readers.
 */
class HypergraphReader {
public:
    virtual ~HypergraphReader() = default;

    /**
     * Given the path of an input file, read the hypergraph contained.
     *
     * @param[in]  path  input file path
     */
    virtual Hypergraph read(std::string_view path) = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_HYPERGRAPH_READER_HPP_
