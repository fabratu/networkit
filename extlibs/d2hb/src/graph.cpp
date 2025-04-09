#include <d2hb/graph.h>

#include <algorithm>

namespace d2hb {
Weight invalidWeight() { return std::numeric_limits<float>::infinity(); }
EdgeID invalidEdgeID() { return std::numeric_limits<EdgeID>::max(); }

Edge invalidEdge() { return Edge(invalidVertex(), invalidTarget()); }

Target invalidTarget() { return {invalidVertex(), {invalidWeight(), invalidEdgeID()}}; }

std::vector<Degree> degrees_from(Edges const& edges) {
    unsigned int count = graph::vertex_count(edges);
    std::vector<Degree> associated_degree(count, 0u);

    for (Edge const& edge : edges) {
        associated_degree[edge.source]++;
    }

    return associated_degree;
}

namespace graph {

Vertex vertex_count(Edges const& edges) {
    // Determine n as the maximal node ID.
    Vertex n = 0;
    for (Edge const& edge : edges) {
        Vertex const vertex = std::max(edge.source, edge.target.vertex);
        n = std::max(n, vertex);
    }

    n += 1;
    return n;
}

Vertex non_zero_vertex_count(Edges const& edges) {
    std::vector<Vertex> vertices;

    for (Edge const& edge : edges) {
        vertices.push_back(edge.source);
        vertices.push_back(edge.target.vertex);
    }

    std::sort(std::begin(vertices), std::end(vertices));
    auto last = std::unique(std::begin(vertices), std::end(vertices));
    vertices.erase(last, std::end(vertices));

    return vertices.size();
}

} // namespace graph

} // namespace d2hb
