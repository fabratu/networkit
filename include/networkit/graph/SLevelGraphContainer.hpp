#ifndef NETWORKIT_GRAPH_S_LEVEL_GRAPH_CONTAINER_HPP_
#define NETWORKIT_GRAPH_S_LEVEL_GRAPH_CONTAINER_HPP_

#include <memory>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/** Builds the cumulative s-level intersection graphs of a hypergraph. */
class SLevelGraphContainer final {
public:
    explicit SLevelGraphContainer(const Hypergraph &hGraph) : hGraph(hGraph) {}

    void build();

    /** Returns one more than the greatest nonempty intersection level. */
    count getMaxLevel() const noexcept { return levels.empty() ? 1 : levels.size() - 1; }

    const Graph &graph(count s) const;
    const std::vector<node> &activeNodes(count s) const;
    count maximumDegree(count s) const;

private:
    struct Level {
        Graph graph;
        std::vector<node> activeNodes;
        count maximumDegree;
    };

    const Level &level(count s) const;

    const Hypergraph &hGraph;
    std::vector<std::shared_ptr<Level>> levels;
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_S_LEVEL_GRAPH_CONTAINER_HPP_
