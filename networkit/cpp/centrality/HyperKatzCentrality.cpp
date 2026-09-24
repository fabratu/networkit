#include <algorithm>
#include <cfloat>
#include <cmath>
#include <stdexcept>

#include <networkit/centrality/HyperKatzCentrality.hpp>

namespace NetworKit {

HyperKatzCentrality::HyperKatzCentrality(const Hypergraph &hGraph, count k, bool groupOnly,
                                         double tolerance)
    : hGraph{hGraph}, matrices{hGraph}, k{k}, groupOnly{groupOnly}, rankTolerance{tolerance} {
    if (k == 0 || k > hGraph.numberOfEdges())
        throw std::invalid_argument("k must be between one and the number of hyperedges");
    if (tolerance < 0)
        throw std::invalid_argument("The ranking tolerance must be non-negative");

    matrices.build(SMatrixType::Level);

    Vector ones(hGraph.upperEdgeIdBound(), 1.0);
    alphaByLevel.resize(matrices.getMaxLevel() + 1, 0.0);
    matrices.forLevels([&](count s, const CSRMatrix &matrix) {
        const Vector degrees = matrix * ones;
        count maxDegree = 0;
        hGraph.forEdges(
            [&](edgeid eid) { maxDegree = std::max(maxDegree, static_cast<count>(degrees[eid])); });

        const double alpha = 1.0 / (static_cast<double>(maxDegree) + 1.0);
        levelMatrices.push_back(&matrix);
        levelMaxDegrees.push_back(maxDegree);
        levelAlphas.push_back(alpha);
        alphaByLevel[s] = alpha;
    });

    if (levelMatrices.empty())
        throw std::runtime_error(
            "At least one non-empty s-level matrix is required for HyperKatzCentrality");
}

void HyperKatzCentrality::run() {
    const count dimension = hGraph.upperEdgeIdBound();

    nPaths.clear();
    nPaths.emplace_back(dimension, 0.0);
    hGraph.forEdges([&](edgeid eid) { nPaths[0][eid] = 1.0; });

    currentPaths.clear();
    currentPaths.reserve(levelMatrices.size());
    for (index i = 0; i < levelMatrices.size(); ++i) {
        currentPaths.emplace_back(dimension, 0.0);
        hGraph.forEdges([&](edgeid eid) { currentPaths.back()[eid] = 1.0; });
    }

    activeRanking.clear();
    activeRanking.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { activeRanking.push_back(eid); });

    scoreData.assign(dimension, 0.0);
    baseData.assign(dimension, 0.0);
    boundData.assign(dimension, DBL_MAX);
    iterationReached = 0;

    do {
        doIteration();
    } while (!checkConvergence());

    hasRun = true;
}

const std::vector<double> &HyperKatzCentrality::scores() const {
    assureFinished();
    return scoreData;
}

double HyperKatzCentrality::score(edgeid eid) const {
    assureFinished();
    return scoreData.at(eid);
}

std::vector<std::pair<edgeid, double>> HyperKatzCentrality::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, scoreData[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

edgeid HyperKatzCentrality::top(count n) const {
    assureFinished();
    return activeRanking.at(n);
}

double HyperKatzCentrality::bound(edgeid eid) const {
    assureFinished();
    return boundData.at(eid);
}

double HyperKatzCentrality::getAlpha(count s) const {
    if (s >= alphaByLevel.size())
        throw std::out_of_range("The s-level exceeds the maximum level");
    if (s == 0 || alphaByLevel[s] == 0.0)
        throw std::invalid_argument("The s-level matrix is empty");
    return alphaByLevel[s];
}

bool HyperKatzCentrality::areDistinguished(edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (scoreData[eid1] < scoreData[eid2])
        std::swap(eid1, eid2);
    return scoreData[eid1] > boundData[eid2];
}

bool HyperKatzCentrality::areSufficientlyRanked(edgeid high, edgeid low) const {
    return scoreData[high] > boundData[low] - rankTolerance;
}

void HyperKatzCentrality::doIteration() {
    const count r = iterationReached + 1;
    const count dimension = hGraph.upperEdgeIdBound();
    nPaths.emplace_back(dimension, 0.0);
    std::vector<double> lowerCorrection(dimension, 0.0);
    std::vector<double> upperCorrection(dimension, 0.0);

    for (index i = 0; i < levelMatrices.size(); ++i) {
        currentPaths[i] = *levelMatrices[i] * currentPaths[i];

        const double alpha = levelAlphas[i];
        const count maxDegree = levelMaxDegrees[i];
        const double alphaPower = std::pow(alpha, static_cast<double>(r));
        const double nextAlphaPower = alpha * alphaPower;
        const double boundFactor = nextAlphaPower * maxDegree / (1.0 - alpha * maxDegree);

        hGraph.forEdges([&](edgeid eid) {
            const double paths = currentPaths[i][eid];
            nPaths[r][eid] += paths;
            baseData[eid] += alphaPower * paths;
            lowerCorrection[eid] += nextAlphaPower * paths;
            upperCorrection[eid] += boundFactor * paths;
        });
    }

    hGraph.parallelForEdges([&](edgeid eid) {
        scoreData[eid] = baseData[eid]; // + lowerCorrection[eid];
        boundData[eid] = baseData[eid] + upperCorrection[eid];
    });

    ++iterationReached;
}

bool HyperKatzCentrality::checkConvergence() {
    if (activeRanking.size() > k) {
        std::partial_sort(
            activeRanking.begin(), activeRanking.begin() + k, activeRanking.end(),
            [&](edgeid eid1, edgeid eid2) { return scoreData[eid1] > scoreData[eid2]; });

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
              [&](edgeid eid1, edgeid eid2) { return scoreData[eid1] > scoreData[eid2]; });
    for (index i = 1; i < activeRanking.size(); ++i) {
        if (!areSufficientlyRanked(activeRanking[i - 1], activeRanking[i]))
            return false;
    }

    return true;
}

} // namespace NetworKit
