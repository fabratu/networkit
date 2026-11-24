#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>

#include <networkit/centrality/HyperKatzCentralityNaiveSumTight.hpp>

namespace NetworKit {

HyperKatzCentralityNaiveSumTight::HyperKatzCentralityNaiveSumTight(const Hypergraph &hGraph,
                                                                   count k, bool groupOnly,
                                                                   double tolerance)
    : hGraph{hGraph}, matrices{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    auto start = std::chrono::high_resolution_clock::now();
    matrices.build(SMatrixType::Level, nullptr, nullptr, true);

    levelIndexByNumber.resize(matrices.getMaxLevel() + 1, std::numeric_limits<index>::max());
    matrices.forLevels([&](count s, const DCSRMatrix &matrix) {
        if (matrix.nnz() == 0)
            return;
        levelIndexByNumber[s] = levelMatrices.size();
        levelMatrices.push_back(&matrix);
        levelAlphas.push_back(matrix.alphaVector());
        levelUpperFactors.push_back(matrix.upperBoundVector(levelAlphas.back()));
    });
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Setup: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms"
              << std::endl;

    if (levelMatrices.empty())
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");
}

void HyperKatzCentralityNaiveSumTight::run() {
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

const Vector &HyperKatzCentralityNaiveSumTight::scores() const {
    assureFinished();
    return msLowerBound;
}

double HyperKatzCentralityNaiveSumTight::score(edgeid eid) const {
    assureFinished();
    return msLowerBound[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentralityNaiveSumTight::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msLowerBound[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentralityNaiveSumTight::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentralityNaiveSumTight::bound(edgeid eid) const {
    assureFinished();
    return msUpperBound[eid];
}

double HyperKatzCentralityNaiveSumTight::getAlpha(count s, edgeid eid) const {
    if (s >= levelIndexByNumber.size())
        throw std::out_of_range("The s-level exceeds the maximum level");
    if (s == 0 || levelIndexByNumber[s] == std::numeric_limits<index>::max())
        throw std::invalid_argument("The s-level matrix is empty");
    if (eid >= hGraph.upperEdgeIdBound() || !hGraph.hasEdge(eid))
        throw std::out_of_range("The hyperedge does not exist");
    return levelAlphas[levelIndexByNumber[s]][eid];
}

bool HyperKatzCentralityNaiveSumTight::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (msLowerBound[eid1] < msLowerBound[eid2])
        std::swap(eid1, eid2);
    return msLowerBound[eid1] > msLowerBound[eid2];
}

bool HyperKatzCentralityNaiveSumTight::areSufficientlyRanked(edgeid high, edgeid low) const {
    return msLowerBound[high] > msUpperBound[low] - rankTolerance;
}

void HyperKatzCentralityNaiveSumTight::doIteration() {
    msUpperBound.parallelForElements(
        [&](index element, double &value) { value = msLowerBound[element]; });

    for (index i = 0; i < levelMatrices.size(); ++i) {
        levelMatrices[i]->multiplyUnitValuesInto(currentPaths[i], nextPaths[i]);
        const Vector &alphas = levelAlphas[i];
        const Vector &upperFactors = levelUpperFactors[i];
#pragma omp parallel for
        for (omp_index row = 0; row < static_cast<omp_index>(currentPaths[i].getDimension());
             ++row) {
            const double pathCount = nextPaths[i][row] * alphas[row];
            nextPaths[i][row] = pathCount;
            msLowerBound[row] += pathCount;
            msUpperBound[row] +=
                (1.0 + upperFactors[row] * std::pow(alphas[row], iterationReached + 1)) * pathCount;
        }
        std::swap(currentPaths[i], nextPaths[i]);
    }
    ++iterationReached;
}

bool HyperKatzCentralityNaiveSumTight::checkConvergence() {
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
