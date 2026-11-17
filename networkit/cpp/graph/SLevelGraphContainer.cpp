#include <algorithm>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <tlx/sort/parallel_mergesort.hpp>

#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/graph/SLevelGraphContainer.hpp>

namespace NetworKit {

namespace {

struct EdgePair {
    edgeid first;
    edgeid second;

    bool operator==(const EdgePair &) const = default;
};

} // namespace

void SLevelGraphContainer::build() {
    levels.clear();

    const count dimension = hGraph.upperEdgeIdBound();
    const count nodeBound = hGraph.upperNodeIdBound();
    std::vector<index> pairOffsets(nodeBound + 1, 0);
    for (node u = 0; u < nodeBound; ++u) {
        if (!hGraph.hasNode(u))
            continue;
        const count degree = hGraph.edgesOf(u).size();
        pairOffsets[u + 1] = degree > 1 ? degree * (degree - 1) / 2 : 0;
    }
    std::partial_sum(pairOffsets.begin(), pairOffsets.end(), pairOffsets.begin());

    const auto buildFromPairs = [&](auto pairs, auto firstOf, auto secondOf, auto less) {
        using Pair = typename decltype(pairs)::value_type;
        if (pairs.size() > 1)
            tlx::parallel_mergesort(pairs.begin(), pairs.end(), less,
                                    std::min<count>(pairs.size(), Aux::getMaxNumberOfThreads()));

        count maximumIntersection = 0;
        std::vector<count> bucketSizes(1, 0);
        for (index begin = 0; begin < pairs.size();) {
            index end = begin + 1;
            while (end < pairs.size() && pairs[end] == pairs[begin])
                ++end;
            const count intersection = end - begin;
            if (bucketSizes.size() <= intersection)
                bucketSizes.resize(intersection + 1, 0);
            ++bucketSizes[intersection];
            maximumIntersection = std::max(maximumIntersection, intersection);
            begin = end;
        }

        std::vector<std::vector<Pair>> exactPairs(maximumIntersection + 1);
        for (count s = 1; s <= maximumIntersection; ++s)
            exactPairs[s].reserve(bucketSizes[s]);
        for (index begin = 0; begin < pairs.size();) {
            index end = begin + 1;
            while (end < pairs.size() && pairs[end] == pairs[begin])
                ++end;
            exactPairs[end - begin].push_back(pairs[begin]);
            begin = end;
        }
        std::vector<Pair>{}.swap(pairs);

        levels.resize(maximumIntersection + 2);
        for (count s = maximumIntersection; s > 0; --s) {
            if (exactPairs[s].empty()) {
                levels[s] = levels[s + 1];
                continue;
            }

            std::vector<count> degrees(dimension, 0);
            count edgeCount = 0;
            for (count level = s; level <= maximumIntersection; ++level) {
                for (const Pair &pair : exactPairs[level]) {
                    const edgeid first = firstOf(pair);
                    const edgeid second = secondOf(pair);
                    ++degrees[first];
                    ++degrees[second];
                    ++edgeCount;
                }
            }
            std::vector<index> rowOffsets(dimension + 1, 0);
            for (node u = 0; u < dimension; ++u)
                rowOffsets[u + 1] = rowOffsets[u] + degrees[u];

            std::vector<index> columns(2 * edgeCount);
            auto positions = rowOffsets;
            for (count level = s; level <= maximumIntersection; ++level) {
                for (const Pair &pair : exactPairs[level]) {
                    const edgeid first = firstOf(pair);
                    const edgeid second = secondOf(pair);
                    columns[positions[first]++] = second;
                    columns[positions[second]++] = first;
                }
            }

            Graph graph(dimension, false, false);
#pragma omp parallel for schedule(guided)
            for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                if (degrees[row] == 0)
                    continue;
                graph.preallocateUndirected(row, degrees[row]);
                for (index entry = rowOffsets[row]; entry < rowOffsets[row + 1]; ++entry)
                    graph.addPartialEdge(unsafe, row, columns[entry]);
            }
            graph.setEdgeCount(unsafe, edgeCount);
            std::vector<node> activeNodes;
            activeNodes.reserve(edgeCount >= dimension / 2 ? dimension : 2 * edgeCount);
            count maximumDegree = 0;
            for (node u = 0; u < dimension; ++u) {
                if (degrees[u] != 0)
                    activeNodes.push_back(u);
                maximumDegree = std::max(maximumDegree, degrees[u]);
            }
            levels[s] = std::make_shared<Level>(
                Level{std::move(graph), std::move(activeNodes), maximumDegree});
        }
    };

    if (dimension <= std::numeric_limits<std::uint32_t>::max()) {
        std::vector<std::uint64_t> pairs(pairOffsets.back());
        hGraph.parallelForNodes([&](node u) {
            index position = pairOffsets[u];
            const auto &incident = hGraph.edgesOf(u);
            for (auto first = incident.begin(); first != incident.end(); ++first) {
                for (auto second = std::next(first); second != incident.end(); ++second) {
                    const auto [eid1, eid2] = std::minmax(*first, *second);
                    pairs[position++] =
                        (static_cast<std::uint64_t>(eid1) << 32) | static_cast<std::uint32_t>(eid2);
                }
            }
        });
        buildFromPairs(
            std::move(pairs), [](std::uint64_t pair) { return static_cast<edgeid>(pair >> 32); },
            [](std::uint64_t pair) {
                return static_cast<edgeid>(static_cast<std::uint32_t>(pair));
            },
            std::less<std::uint64_t>{});
    } else {
        std::vector<EdgePair> pairs(pairOffsets.back());
        hGraph.parallelForNodes([&](node u) {
            index position = pairOffsets[u];
            const auto &incident = hGraph.edgesOf(u);
            for (auto first = incident.begin(); first != incident.end(); ++first) {
                for (auto second = std::next(first); second != incident.end(); ++second) {
                    const auto [eid1, eid2] = std::minmax(*first, *second);
                    pairs[position++] = {eid1, eid2};
                }
            }
        });
        const auto less = [](const EdgePair &lhs, const EdgePair &rhs) {
            return lhs.first < rhs.first || (lhs.first == rhs.first && lhs.second < rhs.second);
        };
        buildFromPairs(
            std::move(pairs), [](const EdgePair &pair) { return pair.first; },
            [](const EdgePair &pair) { return pair.second; }, less);
    }
}

const SLevelGraphContainer::Level &SLevelGraphContainer::level(count s) const {
    if (s == 0 || s >= getMaxLevel())
        throw std::out_of_range("The s-level exceeds the graph container's levels");
    if (levels[s] == nullptr)
        throw std::logic_error("S-level graphs have not been built");
    return *levels[s];
}

const Graph &SLevelGraphContainer::graph(count s) const {
    return level(s).graph;
}

const std::vector<node> &SLevelGraphContainer::activeNodes(count s) const {
    return level(s).activeNodes;
}

count SLevelGraphContainer::maximumDegree(count s) const {
    return level(s).maximumDegree;
}

} // namespace NetworKit
