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

    std::vector<double> levelAlphas;
    matrices.build(SMatrixType::Delta, &currentPaths, &levelAlphas);
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
}

void HyperKatzCentralityDeltaInPlace::run() {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const count dimension = hGraph.upperEdgeIdBound();

    // currentPaths.reset();

    activeRanking.clear();
    activeRanking.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { activeRanking.push_back(eid); });

    msLowerBound = Vector(dimension, 0.0);
    // TODO: reset msUpperCorrection

    lowerCorrection = Vector(matrices.getMaxLevel(), 1.0);

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
    nextPaths.reset();

    matrices.forLevels([&](count s, const ACSRMatrix &matrix) {
        matrix.multiplyInto(currentPaths, nextPaths);
        lowerCorrection[s - 1] *= alphas[s - 1];
        upperCorrection[s - 1] *= alphas[s - 1];
    });

    std::swap(currentPaths, nextPaths);

    // alphas.parallelForElements([&](const int &i, double &element) {
    //     lowerCorrection[i] *= element;
    //     upperCorrection[i] *= element;
    // });

    // double alpha = levelAlphas[s - 1];
    // double alphaPower = std::pow(alpha, static_cast<double>(r));
    // auto deg = levelMaxDegrees[s - 1];
    // TODO: can we do this better with double saving alpha + alphaPower and do vectorized
    // *= and = seperatly?
    //     lowerCorrection[s - 1] *= alpha;
    //     upperCorrection[s - 1] = alphaPower * alpha * deg / (1.0 - alpha * deg);
    // });

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

    // auto lowerSetter = [&](int i, double &element) { element *= levelAlphas[i]; };
    // auto upperSetter = [&](int i, double &element) {

    // };

    // lowerCorrection.forElements(lowerSetter);
    // upperCorrection.parallelForElements(upperSetter);

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
