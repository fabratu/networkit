/*
 * CustomFormatGraphWriter.hpp
 */

// networkit-format

#ifndef NETWORKIT_IO_CUSTOM_FORMAT_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_CUSTOM_FORMAT_GRAPH_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Reader for the new custom file format.
 */
class CustomFormatGraphWriter final : public GraphWriter {

public:
    CustomFormatGraphWriter() = default;

    void write(const Graph &G, const std::string &path) override;

    void utilFunction(const std::string &option);
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_CUSTOM_FORMAT_GRAPH_WRITER_HPP_
