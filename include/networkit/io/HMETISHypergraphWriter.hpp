/*
 * HMETISHypergraphWriter.hpp
 */

#ifndef NETWORKIT_IO_HMETIS_HYPERGRAPH_WRITER_HPP_
#define NETWORKIT_IO_HMETIS_HYPERGRAPH_WRITER_HPP_

#include <string_view>

#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Writer for the hMETIS hypergraph format, including augmented incidence weights.
 */
class HMETISHypergraphWriter final {
public:
    HMETISHypergraphWriter() = default;

    /**
     * Write @a hypergraph to @a path. The format flag is selected from 0, 1, 20, and 21 based on
     * the hypergraph's edge- and incidence-weighted properties.
     *
     * @param hypergraph Hypergraph to write.
     * @param path Output file path.
     */
    void write(const Hypergraph &hypergraph, std::string_view path) const;
};

} // namespace NetworKit

#endif // NETWORKIT_IO_HMETIS_HYPERGRAPH_WRITER_HPP_
