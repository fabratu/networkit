#include <algorithm>
#include <stdexcept>

#include <networkit/centrality/HyperKatzCentralityGraph.hpp>

namespace NetworKit {

HyperKatzCentralityGraph::HyperKatzCentralityGraph(const Hypergraph &hGraph, count k,
                                                   bool groupOnly, double tolerance)
    : hGraph{hGraph}, graphs{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    graphs.build();
    alphaByLevel.resize(graphs.getMaxLevel() + 1, 0.0);
    for (count s = 1; s < graphs.getMaxLevel(); ++s) {
        const count maxDegree = graphs.maximumDegree(s);
        if (maxDegree == 0)
            continue;
        const double alpha = 1.0 / (static_cast<double>(maxDegree) + 1.0);
        levelNumbers.push_back(s);
        levelAlphas.push_back(alpha);
        levelUpperFactors.push_back(alpha * maxDegree / (1.0 - alpha * maxDegree));
        alphaByLevel[s] = alpha;
    }

    if (levelNumbers.empty())
        throw std::runtime_error(
            "At least one non-empty s-level graph is required for HyperKatzCentrality");
}

void HyperKatzCentralityGraph::run() {
    const count dimension = hGraph.upperEdgeIdBound();
    currentPaths.clear();
    nextPaths.clear();
    currentPaths.reserve(levelNumbers.size());
    nextPaths.reserve(levelNumbers.size());
    for (count s : levelNumbers) {
        currentPaths.emplace_back(dimension, 0.0);
        nextPaths.emplace_back(dimension, 0.0);
        for (node u : graphs.activeNodes(s))
            currentPaths.back()[u] = 1.0;
    }
    alphaPowers = levelAlphas;

    activeRanking.clear();
    activeRanking.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { activeRanking.push_back(eid); });

    msLowerBound = Vector(dimension, 0.0);
    msUpperBound = Vector(dimension, 0.0);
    iterationReached = 0;

    do {
        doIteration();
    } while (!checkConvergence());
    hasRun = true;
}

const Vector &HyperKatzCentralityGraph::scores() const {
    assureFinished();
    return msLowerBound;
}

double HyperKatzCentralityGraph::score(edgeid eid) const {
    assureFinished();
    return msLowerBound[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentralityGraph::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msLowerBound[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentralityGraph::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentralityGraph::bound(edgeid eid) const {
    assureFinished();
    return msUpperBound[eid];
}

bool HyperKatzCentralityGraph::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (msLowerBound[eid1] < msLowerBound[eid2])
        std::swap(eid1, eid2);
    return msLowerBound[eid1] > msLowerBound[eid2];
}

double HyperKatzCentralityGraph::getAlpha(count s) const {
    if (s >= alphaByLevel.size())
        throw std::out_of_range("The s-level exceeds the maximum level");
    if (s == 0 || alphaByLevel[s] == 0.0)
        throw std::invalid_argument("The s-level graph is empty");
    return alphaByLevel[s];
}

bool HyperKatzCentralityGraph::areSufficientlyRanked(edgeid high, edgeid low) const {
    return msLowerBound[high] > msUpperBound[low] - rankTolerance;
}

void HyperKatzCentralityGraph::doIteration() {
    msUpperBound.parallelForElements(
        [&](index element, double &value) { value = msLowerBound[element]; });

    for (index i = 0; i < levelNumbers.size(); ++i) {
        const Graph &graph = graphs.graph(levelNumbers[i]);
        const auto &active = graphs.activeNodes(levelNumbers[i]);
        const double lowerFactor = alphaPowers[i];
        const double upperFactor = lowerFactor * levelUpperFactors[i];
        const auto &current = currentPaths[i];
        auto &next = nextPaths[i];

#pragma omp parallel for schedule(guided)
        for (omp_index position = 0; position < static_cast<omp_index>(active.size()); ++position) {
            const node u = active[position];
            double paths = 0.0;
            graph.forNeighborsOf(u, [&](node v) { paths += current[v]; });
            next[u] = paths;

            const double lowerContribution = lowerFactor * paths;
            msLowerBound[u] += lowerContribution;
            msUpperBound[u] += lowerContribution + upperFactor * paths;
        }
        std::swap(currentPaths[i], nextPaths[i]);
        alphaPowers[i] *= levelAlphas[i];
    }
    ++iterationReached;
}

bool HyperKatzCentralityGraph::checkConvergence() {
    if (activeRanking.size() > k) {
        std::partial_sort(
            activeRanking.begin(), activeRanking.begin() + k, activeRanking.end(),
            [&](edgeid eid1, edgeid eid2) { return msLowerBound[eid1] > msLowerBound[eid2]; });

        const edgeid kth = activeRanking[k - 1];
        activeRanking.erase(
            std::remove_if(activeRanking.begin() + k, activeRanking.end(),
                           [&](edgeid eid) { return areSufficientlyRanked(kth, eid); }),
            activeRanking.end());
    }

    if (activeRanking.size() > k)
        return false;
    if (groupOnly)
        return true;

    std::sort(activeRanking.begin(), activeRanking.end(),
              [&](edgeid eid1, edgeid eid2) { return msLowerBound[eid1] > msLowerBound[eid2]; });
    for (index j = 1; j < activeRanking.size(); ++j) {
        if (!areSufficientlyRanked(activeRanking[j - 1], activeRanking[j]))
            return false;
    }
    return true;
}

} // namespace NetworKit
