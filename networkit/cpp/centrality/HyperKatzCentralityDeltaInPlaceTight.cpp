#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include <networkit/centrality/HyperKatzCentralityDeltaInPlaceTight.hpp>

namespace NetworKit {

HyperKatzCentralityDeltaInPlaceTight::HyperKatzCentralityDeltaInPlaceTight(const Hypergraph &hGraph, count k,
                                                                 bool groupOnly, double tolerance)
    : hGraph{hGraph}, matrices{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    auto start = std::chrono::high_resolution_clock::now();
    matrices.build(SMatrixType::Delta, &currentPaths, nullptr, true);
    nextPaths = currentPaths;
    nextPaths.reset();
    if (matrices.getMaxLevel() <= 1 || matrices.getMatrix(1).nnz() == 0)
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");

    levelAlphas = currentPaths;
    levelUpperFactors = currentPaths;
    std::vector<count> degrees(hGraph.upperEdgeIdBound(), 0);
    for (count s = matrices.getMaxLevel() - 1; s > 0; --s) {
        matrices.getMatrix(s).parallelForRowDegrees(
            [&](node u, count exactDegree) { degrees[u] += exactDegree; });
#pragma omp parallel for schedule(static)
        for (omp_index row = 0; row < static_cast<omp_index>(hGraph.numberOfEdges()); ++row) {
            if (currentPaths.rowLength(row) < s)
                continue;
            const double degree = static_cast<double>(degrees[row]);
            const double alpha = 1.0 / (degree + 1.0);
            levelAlphas.rowData(row)[s - 1] = alpha;
            levelUpperFactors.rowData(row)[s - 1] =
                degree * alpha / (1.0 - alpha * degree);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Setup: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
}

void HyperKatzCentralityDeltaInPlaceTight::run() {
    auto start = std::chrono::high_resolution_clock::now();
    const count dimension = hGraph.upperEdgeIdBound();

#pragma omp parallel for schedule(static)
    for (omp_index row = 0; row < static_cast<omp_index>(hGraph.numberOfEdges()); ++row)
        std::fill_n(currentPaths.rowData(row), currentPaths.rowLength(row), 1.0);
    nextPaths.reset();

    activeRanking.clear();
    activeRanking.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { activeRanking.push_back(eid); });

    msLowerBound = Vector(dimension, 0.0);
    msUpperBound = Vector(dimension, 0.0);

    iterationReached = 1;

    do {
        doIteration();
    } while (!checkConvergence());

    hasRun = true;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Run: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;
}

const Vector &HyperKatzCentralityDeltaInPlaceTight::scores() const {
    assureFinished();
    return msLowerBound;
}

double HyperKatzCentralityDeltaInPlaceTight::score(edgeid eid) const {
    assureFinished();
    return msLowerBound[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentralityDeltaInPlaceTight::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msLowerBound[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentralityDeltaInPlaceTight::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentralityDeltaInPlaceTight::bound(edgeid eid) const {
    assureFinished();
    return msUpperBound[eid];
}

bool HyperKatzCentralityDeltaInPlaceTight::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (msLowerBound[eid1] < msLowerBound[eid2])
        std::swap(eid1, eid2);
    return msLowerBound[eid1] > msLowerBound[eid2];
}

bool HyperKatzCentralityDeltaInPlaceTight::areSufficientlyRanked(edgeid high, edgeid low) const {
    return msLowerBound[high] > msUpperBound[low] - rankTolerance;
}

void HyperKatzCentralityDeltaInPlaceTight::doIteration() {
    nextPaths.reset();
    matrices.multiplyInto(currentPaths, nextPaths);
    msUpperBound.parallelForElements(
        [&](index element, double &value) { value = msLowerBound[element]; });

    const count iteration = iterationReached + 1;
#pragma omp parallel for schedule(static)
    for (omp_index row = 0; row < static_cast<omp_index>(hGraph.numberOfEdges()); ++row) {
        const double *const alphas = levelAlphas.rowData(row);
        const double *const upperFactors = levelUpperFactors.rowData(row);
        double *const paths = nextPaths.rowData(row);
        double lower = msLowerBound[row];
        double upper = msUpperBound[row];
        for (count level = 0; level < nextPaths.rowLength(row); ++level) {
            const double pathCount = paths[level] * alphas[level];
            paths[level] = pathCount;
            lower += pathCount;
            upper +=
                (1.0 + upperFactors[level] * std::pow(alphas[level], iteration)) * pathCount;
        }
        msLowerBound[row] = lower;
        msUpperBound[row] = upper;
    }

    std::swap(currentPaths, nextPaths);

    ++iterationReached;
}

bool HyperKatzCentralityDeltaInPlaceTight::checkConvergence() {
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
