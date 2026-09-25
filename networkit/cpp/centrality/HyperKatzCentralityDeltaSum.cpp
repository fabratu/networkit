#include <algorithm>
#include <cfloat>
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

    matrices.build(SMatrixType::Delta);

    Vector ones(hGraph.upperEdgeIdBound(), 1.0);
    // double levelTol = rankTolerance / matrices.getMaxLevel();
    matrices.forLevels([&](count s, const CSRMatrix &matrix) {
        const Vector degrees = matrix * ones;
        count maxDegree = 0;
        hGraph.forEdges(
            [&](edgeid eid) { maxDegree = std::max(maxDegree, static_cast<count>(degrees[eid])); });

        deltaMatrices.push_back(&matrix);
        levelMaxDegrees.push_back(maxDegree);
        // Add degree of delta matrix to all lower levels
        for (int i = 0; i < levelMaxDegrees.size() - 1; i++) {
            levelMaxDegrees[i] += maxDegree;
        }
    });

    for (auto val : levelMaxDegrees) {
        const double alpha = 1.0 / (static_cast<double>(val) + 1.0);
        levelAlphas.push_back(alpha);
    }

    if (deltaMatrices.empty())
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");
}

void HyperKatzCentralityDeltaSum::run() {
    const count dimension = hGraph.upperEdgeIdBound();

    currentPaths.clear();
    currentPaths.reserve(deltaMatrices.size());
    for (index i = 0; i < deltaMatrices.size(); ++i) {
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

    for (index i = 0; i < deltaMatrices.size(); ++i) {
        deltaCorrection.fill(0.0);
        for (index j = 0; j <= i; j++) {
            deltaCorrection += *deltaMatrices[j] * currentPaths[i];
        }
        currentPaths[i] = deltaCorrection;
        const double alpha = levelAlphas[i];
        const count maxDegree = levelMaxDegrees[i];
        const double alphaPower = std::pow(alpha, static_cast<double>(r));
        const double nextAlphaPower = alpha * alphaPower;
        const double boundFactor = maxDegree / (1.0 - alpha * maxDegree);
        msLowerBound += alphaPower * currentPaths[i];
        msUpperBound += nextAlphaPower * boundFactor * currentPaths[i];
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
