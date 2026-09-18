/*
 * HMETISHypergraphReader.hpp
 */

#ifndef NETWORKIT_IO_HMETIS_HYPERGRAPH_READER_HPP_
#define NETWORKIT_IO_HMETIS_HYPERGRAPH_READER_HPP_

#include <string_view>

#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the hMETIS hypergraph format.
 *
 * Besides the standard unweighted (0) and hyperedge-weighted (1) variants, this reader supports
 * the augmented incidence-weighted variants 20 and 21. In these variants, every node id is
 * followed by its incidence weight. Node ids in the file are one-based.
 */
class HMETISHypergraphReader final {
public:
    HMETISHypergraphReader() = default;

    /**
     * Read a hypergraph from @a path.
     *
     * @param path Input file path.
     * @return The hypergraph read from the file.
     */
    Hypergraph read(std::string_view path) const;
};

} // namespace NetworKit

#endif // NETWORKIT_IO_HMETIS_HYPERGRAPH_READER_HPP_
