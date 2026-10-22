#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/HyperKatzCentralityDeltaInPlace.hpp>

namespace NetworKit {

HyperKatzCentralityDeltaInPlace::HyperKatzCentralityDeltaInPlace(const Hypergraph &hGraph, count k,
                                                                 bool groupOnly, double tolerance)
    : hGraph{hGraph}, matrices{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> levelAlphas;
    matrices.build(SMatrixType::Delta, &currentPaths, &levelAlphas, true);
    nextPaths = currentPaths;
    nextPaths.reset();

    alphas = Vector(levelAlphas);
    upperCorrection = Vector(levelAlphas.size(), 0.0);
    for (index i = 0; i < levelAlphas.size(); ++i) {
        const double alpha = levelAlphas[i];
        const double maxDegree = std::round(1.0 / alpha - 1.0);
        upperCorrection[i] = alpha * maxDegree / (1.0 - alpha * maxDegree);
    }

    if (matrices.numberOfMatrices() == 0)
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Setup: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
}

void HyperKatzCentralityDeltaInPlace::run() {
    auto start = std::chrono::high_resolution_clock::now();
    const count dimension = hGraph.upperEdgeIdBound();

    activeRanking.clear();
    activeRanking.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { activeRanking.push_back(eid); });

    msLowerBound = Vector(dimension, 0.0);

    lowerCorrection = Vector(matrices.getMaxLevel(), 1.0);

    iterationReached = 0;

    do {
        doIteration();
    } while (!checkConvergence());

    hasRun = true;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Run: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
}

const Vector &HyperKatzCentralityDeltaInPlace::scores() const {
    assureFinished();
    return msLowerBound;
}

double HyperKatzCentralityDeltaInPlace::score(edgeid eid) const {
    assureFinished();
    return msLowerBound[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentralityDeltaInPlace::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msLowerBound[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentralityDeltaInPlace::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentralityDeltaInPlace::bound(edgeid eid) const {
    assureFinished();
    return msUpperBound[eid];
}

bool HyperKatzCentralityDeltaInPlace::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (msLowerBound[eid1] < msLowerBound[eid2])
        std::swap(eid1, eid2);
    return msLowerBound[eid1] > msLowerBound[eid2];
}

bool HyperKatzCentralityDeltaInPlace::areSufficientlyRanked(edgeid high, edgeid low) const {
    return msLowerBound[high] > msUpperBound[low] - rankTolerance;
}

void HyperKatzCentralityDeltaInPlace::doIteration() {
    nextPaths.reset();

    matrices.multiplyInto(currentPaths, nextPaths);
    for (index level = 0; level < matrices.getMaxLevel(); ++level) {
        lowerCorrection[level] *= alphas[level];
        upperCorrection[level] *= alphas[level];
    }

    std::swap(currentPaths, nextPaths);

    msLowerBound += currentPaths * lowerCorrection;
    msUpperBound = currentPaths * upperCorrection;
    msUpperBound += msLowerBound;

    ++iterationReached;
}

bool HyperKatzCentralityDeltaInPlace::checkConvergence() {
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

    if (activeRanking.size() > k) {
        return false;
    }

    if (groupOnly)
        return true;

    std::sort(activeRanking.begin(), activeRanking.end(),
              [&](edgeid eid1, edgeid eid2) { return msLowerBound[eid1] > msLowerBound[eid2]; });
    for (index j = 1; j < activeRanking.size(); ++j) {
        if (!areSufficientlyRanked(activeRanking[j - 1], activeRanking[j])) {
            return false;
        }
    }

    return true;
}

} // namespace NetworKit
