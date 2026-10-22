#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/HyperKatzCentralityDeltaSum.hpp>

namespace NetworKit {

HyperKatzCentralityDeltaSum::HyperKatzCentralityDeltaSum(const Hypergraph &hGraph, count k,
                                                         bool groupOnly, double tolerance)
    : hGraph{hGraph}, matrices{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    auto start = std::chrono::high_resolution_clock::now();
    matrices.build(SMatrixType::Delta, nullptr, &levelAlphas, true);

    bool hasInteractions = false;
    matrices.forLevels(
        [&](count, const DCSRMatrix &matrix) { hasInteractions |= matrix.nnz() != 0; });
    if (!hasInteractions)
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");

    levelMaxDegrees.reserve(levelAlphas.size());
    for (double alpha : levelAlphas)
        levelMaxDegrees.push_back(static_cast<count>(std::round(1.0 / alpha - 1.0)));

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Setup: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
}

void HyperKatzCentralityDeltaSum::run() {
    auto start = std::chrono::high_resolution_clock::now();
    const count dimension = hGraph.upperEdgeIdBound();

    currentPaths.clear();
    currentPaths.reserve(matrices.getMaxLevel());
    for (index i = 0; i < matrices.getMaxLevel(); ++i) {
        currentPaths.emplace_back(dimension, 1.0);
    }

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
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Run: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
}

const Vector &HyperKatzCentralityDeltaSum::scores() const {
    assureFinished();
    return msLowerBound;
}

double HyperKatzCentralityDeltaSum::score(edgeid eid) const {
    assureFinished();
    return msLowerBound[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentralityDeltaSum::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msLowerBound[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentralityDeltaSum::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentralityDeltaSum::bound(edgeid eid) const {
    assureFinished();
    return msUpperBound[eid];
}

bool HyperKatzCentralityDeltaSum::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (msLowerBound[eid1] < msLowerBound[eid2])
        std::swap(eid1, eid2);
    return msLowerBound[eid1] > msLowerBound[eid2];
}

bool HyperKatzCentralityDeltaSum::areSufficientlyRanked(edgeid high, edgeid low) const {
    return msLowerBound[high] > msUpperBound[low] - rankTolerance;
}

void HyperKatzCentralityDeltaSum::doIteration() {
    const count r = iterationReached + 1;
    const count dimension = hGraph.upperEdgeIdBound();
    msUpperBound = Vector(dimension);
    Vector deltaCorrection = Vector(dimension);

    for (index level = 1; level <= matrices.getMaxLevel(); ++level) {
        deltaCorrection.fill(0.0);
        matrices.forLevelsReverseUntil(level, [&](count, const DCSRMatrix &matrix) {
            if (matrix.nnz() > 0)
                matrix.addUnitValuesInto(currentPaths[level - 1], deltaCorrection);
        });

        currentPaths[level - 1] = deltaCorrection;
        const double alpha = levelAlphas[level - 1];
        const count maxDegree = levelMaxDegrees[level - 1];
        const double alphaPower = std::pow(alpha, static_cast<double>(r));
        const double nextAlphaPower = alpha * alphaPower;
        const double boundFactor = maxDegree / (1.0 - alpha * maxDegree);
        msLowerBound += alphaPower * currentPaths[level - 1];
        msUpperBound += nextAlphaPower * boundFactor * currentPaths[level - 1];
    }

    msUpperBound += msLowerBound;

    ++iterationReached;
}

// bool HyperKatzCentralityNaiveSumPerLevelEps::checkGlobalConvergence() {

//     if (activeLevel == 0)
//         return true;

//     for (index i = 0; i < levelMatrices.size(); ++i) {
//         // TODO: make checkconvergence check each level (?), maybe more efficient to check in
//         // doIteration for upperCorrection vs eps_s for every entry
//         checkConvergence;
//     }
// }

bool HyperKatzCentralityDeltaSum::checkConvergence() {
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
