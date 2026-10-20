#include <algorithm>
#include <chrono>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/HyperKatzCentralityNaiveSum.hpp>

namespace NetworKit {

HyperKatzCentralityNaiveSum::HyperKatzCentralityNaiveSum(const Hypergraph &hGraph, count k,
                                                         bool groupOnly, double tolerance)
    : hGraph{hGraph}, matrices{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    auto start = std::chrono::high_resolution_clock::now();
    matrices.build(SMatrixType::Level, nullptr, nullptr, true);

    alphaByLevel.resize(matrices.getMaxLevel() + 1, 0.0);
    matrices.forLevels([&](count s, const DCSRMatrix &matrix) {
        if (matrix.nnz() == 0)
            return;
        const count maxDegree = matrix.maximumRowDegree();
        const double alpha = 1.0 / (static_cast<double>(maxDegree) + 1.0);
        levelMatrices.push_back(&matrix);
        levelAlphas.push_back(alpha);
        levelUpperFactors.push_back(alpha * maxDegree / (1.0 - alpha * maxDegree));
        alphaByLevel[s] = alpha;
    });
    auto end = std::chrono::high_resolution_clock::now();
    INFO("Setup: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());

    if (levelMatrices.empty())
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");
}

void HyperKatzCentralityNaiveSum::run() {
    auto start = std::chrono::high_resolution_clock::now();
    const count dimension = hGraph.upperEdgeIdBound();

    currentPaths.clear();
    currentPaths.reserve(levelMatrices.size());
    nextPaths.clear();
    nextPaths.reserve(levelMatrices.size());
    for (index i = 0; i < levelMatrices.size(); ++i) {
        currentPaths.emplace_back(dimension, 1.0);
        nextPaths.emplace_back(dimension, 0.0);
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
    auto end = std::chrono::high_resolution_clock::now();
    INFO("Run: ", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
}

const Vector &HyperKatzCentralityNaiveSum::scores() const {
    assureFinished();
    return msLowerBound;
}

double HyperKatzCentralityNaiveSum::score(edgeid eid) const {
    assureFinished();
    return msLowerBound[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentralityNaiveSum::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msLowerBound[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentralityNaiveSum::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentralityNaiveSum::bound(edgeid eid) const {
    assureFinished();
    return msUpperBound[eid];
}

double HyperKatzCentralityNaiveSum::getAlpha(count s) const {
    if (s >= alphaByLevel.size())
        throw std::out_of_range("The s-level exceeds the maximum level");
    if (s == 0 || alphaByLevel[s] == 0.0)
        throw std::invalid_argument("The s-level matrix is empty");
    return alphaByLevel[s];
}

bool HyperKatzCentralityNaiveSum::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (msLowerBound[eid1] < msLowerBound[eid2])
        std::swap(eid1, eid2);
    return msLowerBound[eid1] > msLowerBound[eid2];
}

bool HyperKatzCentralityNaiveSum::areSufficientlyRanked(edgeid high, edgeid low) const {
    return msLowerBound[high] > msUpperBound[low] - rankTolerance;
}

void HyperKatzCentralityNaiveSum::doIteration() {
    msUpperBound.parallelForElements(
        [&](index element, double &value) { value = msLowerBound[element]; });

    for (index i = 0; i < levelMatrices.size(); ++i) {
        levelMatrices[i]->multiplyUnitValuesInto(currentPaths[i], nextPaths[i]);
        std::swap(currentPaths[i], nextPaths[i]);

        const double lowerFactor = alphaPowers[i];
        const double upperFactor = lowerFactor * levelUpperFactors[i];
        currentPaths[i].parallelForElements([&](index element, double pathCount) {
            const double lowerContribution = lowerFactor * pathCount;
            msLowerBound[element] += lowerContribution;
            msUpperBound[element] += lowerContribution + upperFactor * pathCount;
        });
        alphaPowers[i] *= levelAlphas[i];
    }
    ++iterationReached;
}

bool HyperKatzCentralityNaiveSum::checkConvergence() {
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
        if (!areSufficientlyRanked(activeRanking[j - 1], activeRanking[j])) {
            return false;
        }
    }

    return true;
}

} // namespace NetworKit
