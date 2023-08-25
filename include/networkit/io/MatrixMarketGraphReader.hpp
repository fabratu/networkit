#ifndef NETWORKIT_IO_MATRIX_MARKET_GRAPH_READER_HPP_
#define NETWORKIT_IO_MATRIX_MARKET_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reader for the matrix market file format documented in
 * https://networkrepository.com/mtx-matrix-market-format.html
 *
 * Does not allow complex fields.
 *
 * TODO: check possible combinations of field values e.g.:
 *       - if field==pattern then coordinate==format
 *       - if symmetry==skew-symmetric then diagonal entries must not be listed
 */
class MatrixMarketGraphReader final : public GraphReader {
public:
    MatrixMarketGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph.
     *
     * @param path
     * @return Graph
     */
    Graph read(const std::string &path) override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_MATRIX_MARKET_GRAPH_READER_HPP_