/*
 * METISGraphReader.cpp
 */

#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/io/CustomFormatGraphReader.hpp>

namespace NetworKit {

Graph CustomFormatGraphReader::read(const std::string& path) {

    GraphBuilder b(10, false);
    auto G = b.toGraph(false);

    /*
     * Doing stuff here ...
     */

    return G;
}

} /* namespace NetworKit */
