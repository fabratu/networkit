#pragma once

#include <limits>
#include <random>
#include <tuple>
#include <vector>

namespace d2hb {

using Weight = double;
using Vertex = uint64_t;
using Degree = Vertex;
using EdgeID = uint64_t;

struct EdgeData {
    Weight weight;
    EdgeID id;
};

struct Target {
    Vertex vertex;
    EdgeData data;
};

struct Edge {
    Edge() = default;

    Edge(Vertex _source, Target _target) : source(_source), target(_target) {}

    Vertex source;
    Target target;
};

using Edges = std::vector<Edge>;
using Targets = std::vector<Target>;

constexpr Vertex invalidVertex() { return std::numeric_limits<Vertex>::max(); }
Weight invalidWeight();
EdgeID invalidEdgeID();
Edge invalidEdge();
Target invalidTarget();

std::vector<Degree> degrees_from(Edges const& edges);

namespace graph {
Vertex vertex_count(Edges const& edges);
Vertex non_zero_vertex_count(Edges const& edges);
Edge into(Vertex source, Target const&);
Target into(Edge const&);
} // namespace graph
} // namespace d2hb
