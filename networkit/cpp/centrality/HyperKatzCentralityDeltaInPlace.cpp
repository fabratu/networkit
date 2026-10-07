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

    matrices.build(SMatrixType::Delta);

    count maxDegree = 0.0;

    SMatrixContainer<DCSRMatrix> maxDegBuilder(hGraph);
    maxDegBuilder.build(SMatrixType::Level);
    Vector ones(hGraph.upperEdgeIdBound(), 1.0);

    // This actually overestimates the maxDegree. For correct result, we need the level matrices (or
    // store the correct value in SMatrixContainer when building delta matrices)
    matrices.forLevels([&](count s, const ACSRMatrix &matrix) {
        const NetworKit::DCSRMatrix &levelRef = maxDegBuilder.getMatrix(s);
        const Vector degrees = levelRef * ones;
        maxDegree = degrees.max();

        // deltaMatrices.push_back(&matrix);
        levelMaxDegrees.push_back(maxDegree);
        const double alpha = 1.0 / (static_cast<double>(maxDegree) + 1.0);
        levelAlphas.push_back(alpha);
        // Add degree of delta matrix to all lower levels
        // for (int i = 0; i < levelMaxDegrees.size() - 1; i++) {
        //     levelMaxDegrees[i] += maxDegree;
        // }
    });

    if (matrices.numberOfMatrices() == 0)
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");
}

void HyperKatzCentralityDeltaInPlace::run() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const count dimension = hGraph.upperEdgeIdBound();

    currentPaths = DenseMatrix(dimension, matrices.getMaxLevel(), 0.0);

    activeRanking.clear();
    activeRanking.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { activeRanking.push_back(eid); });

    msLowerBound = Vector(dimension, 0.0);
    // msUpperBound = Vector(dimension, 0.0);

    lowerCorrection = Vector(matrices.getMaxLevel(), 0.0);
    upperCorrection = Vector(matrices.getMaxLevel(), 0.0);
    auto lowerSetter = [&](int i, double &element) { element = levelAlphas[i]; };
    auto upperSetter = [&](int i, double &element) {
        double alpha = levelAlphas[i];
        double deg = levelMaxDegrees[i];
        element = std::pow(alpha, 2.0) * deg / (1.0 - alpha * deg);
    };

    lowerCorrection.parallelForElements(lowerSetter);
    upperCorrection.parallelForElements(upperSetter);

    iterationReached = 0;

    do {
        doIteration();
    } while (!checkConvergence());

    hasRun = true;
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
    const count r = iterationReached + 1;
    const count dimension = hGraph.upperEdgeIdBound();
    matrices.getMatrix(matrices.getMaxLevel()).resetOther(currentPaths);

    matrices.forLevelsMutable(
        [&](count s, ACSRMatrix &matrix) { matrix.updateOther(currentPaths); });

    msLowerBound += currentPaths * lowerCorrection;
    msUpperBound = currentPaths * upperCorrection;
    msUpperBound += msLowerBound;

    // INFO("LowerBound: ", msLowerBound);
    // INFO("UpperBound: ", msUpperBound);
    // INFO("Iter: ", r);
    // INFO("LowerCorrection: ", lowerCorrection);
    // INFO("UpperCorrection: ", upperCorrection);
    // INFO("Current paths first edge: ", currentPaths.row(0));
    // INFO("Current paths last edge: ", currentPaths.row(currentPaths.numberOfRows() - 1));

    auto lowerSetter = [&](int i, double &element) { element *= levelAlphas[i]; };
    auto upperSetter = [&](int i, double &element) {
        double alpha = levelAlphas[i];
        double alphaPower = std::pow(alpha, static_cast<double>(r + 1));
        double deg = levelMaxDegrees[i];
        element = alphaPower * alpha * deg / (1.0 - alpha * deg);
    };

    lowerCorrection.parallelForElements(lowerSetter);
    upperCorrection.parallelForElements(upperSetter);

    matrices.forLevelsMutable([&](count s, ACSRMatrix &matrix) { matrix.assign(currentPaths); });

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
