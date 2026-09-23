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

    // nPaths.clear();
    // nPaths.emplace_back(dimension, 0.0);

    // hGraph.forEdges([&](edgeid eid) { currentPaths[eid] = 1.0; });

    currentPaths.clear();
    currentPaths.reserve(levelMatrices.size());
    lowerBound.clear();
    lowerBound.reserve(levelMatrices.size());
    upperBound.clear();
    upperBound.reserve(levelMatrices.size());
    for (index i = 0; i < levelMatrices.size(); ++i) {
        currentPaths.emplace_back(dimension, 1.0);
        lowerBound.emplace_back(dimension, 0.0);
        upperBound.emplace_back(dimension, DBL_MAX);
        // hGraph.forEdges([&](edgeid eid) { currentPaths.back()[eid] = 1.0; });
    }

    activeRanking.clear();
    activeRanking.reserve(levelMatrices.size());
    activeLevel.clear();
    activeLevel.reserve(levelMatrices.size());
    for (index i = 0; i < levelMatrices.size(); ++i) {
        activeLevel[i] = true;
        activeRanking[i].reserve(hGraph.numberOfEdges());
        hGraph.forEdges([&](edgeid eid) { activeRanking[i].push_back(eid); });
    }

    msScores = Vector(dimension, 0.0);
    // baseData.assign(dimension, 0.0);
    // boundData.assign(dimension, DBL_MAX);
    iterationReached = 0;

    do {
        doIteration();
    } while (!checkGlobalConvergence());

    hasRun = true;
}

const Vector &HyperKatzCentrality::scores() const {
    assureFinished();
    return msScores;
}

double HyperKatzCentrality::score(edgeid eid) const {
    assureFinished();
    return msScores[eid];
}

std::vector<std::pair<edgeid, double>> HyperKatzCentrality::ranking() const {
    assureFinished();
    std::vector<std::pair<edgeid, double>> result;
    result.reserve(hGraph.numberOfEdges());
    hGraph.forEdges([&](edgeid eid) { result.emplace_back(eid, msScores[eid]); });
    std::sort(result.begin(), result.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first < rhs.first);
    });
    return result;
}

// edgeid HyperKatzCentrality::top(count n) const {
//     assureFinished();
//     return activeRanking.at(n);
// }

// double HyperKatzCentrality::bound(edgeid eid) const {
//     assureFinished();
//     return msBounds[eid];
// }

double HyperKatzCentrality::getAlpha(count s) const {
    if (s >= alphaByLevel.size())
        throw std::out_of_range("The s-level exceeds the maximum level");
    if (s == 0 || alphaByLevel[s] == 0.0)
        throw std::invalid_argument("The s-level matrix is empty");
    return alphaByLevel[s];
}

bool HyperKatzCentrality::areDistinguished(index i, edgeid eid1, edgeid eid2) const {
    assureFinished();
    if (lowerBound[i][eid1] < lowerBound[i][eid2])
        std::swap(eid1, eid2);
    return lowerBound[i][eid1] > lowerBound[i][eid2];
}

bool HyperKatzCentrality::areSufficientlyRanked(index i, edgeid high, edgeid low) const {
    return lowerBound[i][high] > upperBound[i][low] - rankTolerance;
}

// void HyperKatzCentrality::doIteration() {
//     const count r = iterationReached + 1;
//     const count dimension = hGraph.upperEdgeIdBound();
//     nPaths.emplace_back(dimension, 0.0);
//     std::vector<double> lowerCorrection(dimension, 0.0);
//     std::vector<double> upperCorrection(dimension, 0.0);

//     for (index i = 0; i < levelMatrices.size(); ++i) {
//         currentPaths[i] = *levelMatrices[i] * currentPaths[i];

//         const double alpha = levelAlphas[i];
//         const count maxDegree = levelMaxDegrees[i];
//         const double alphaPower = std::pow(alpha, static_cast<double>(r));
//         const double nextAlphaPower = alpha * alphaPower;
//         const double boundFactor = nextAlphaPower * maxDegree / (1.0 - alpha * maxDegree);

//         hGraph.forEdges([&](edgeid eid) {
//             const double paths = currentPaths[i][eid];
//             nPaths[r][eid] += paths;
//             baseData[eid] += alphaPower * paths;
//             lowerCorrection[eid] += nextAlphaPower * paths;
//             upperCorrection[eid] += boundFactor * paths;
//         });
//     }

//     hGraph.parallelForEdges([&](edgeid eid) {
//         scoreData[eid] = baseData[eid] + lowerCorrection[eid];
//         boundData[eid] = baseData[eid] + upperCorrection[eid];
//     });

//     ++iterationReached;
// }

// NOTES:
// - currentPaths holds vectors, currently vector values, but without alpha paths are num paths are
// uints
// - parrallelize over levelMatrices, maybe via parallelForLevel in matrix container

void HyperKatzCentrality::doIteration() {
    const count r = iterationReached + 1;
    const count dimension = hGraph.upperEdgeIdBound();
    // std::vector<double> lowerCorrection(dimension, 0.0);
    // std::vector<double> upperCorrection(dimension, 0.0);

    for (index i = 0; i < levelMatrices.size(); ++i) {
        currentPaths[i] = *levelMatrices[i] * currentPaths[i];

        const double alpha = levelAlphas[i];
        const count maxDegree = levelMaxDegrees[i];
        const double alphaPower = std::pow(alpha, static_cast<double>(r));
        const double nextAlphaPower = alpha * alphaPower;
        const double boundFactor = nextAlphaPower * maxDegree / (1.0 - alpha * maxDegree);

        lowerBound[i] += alphaPower * currentPaths[i];
        upperBound[i] += nextAlphaPower * boundFactor * currentPaths[i];

        //     hGraph.forEdges([&](edgeid eid) {
        //         const double paths = currentPaths[i][eid];
        //         nPaths[r][eid] += paths;
        //         baseData[eid] += alphaPower * paths;
        //         lowerCorrection[eid] += nextAlphaPower * paths;
        //         upperCorrection[eid] += boundFactor * paths;
        //     });
        // }

        // hGraph.parallelForEdges([&](edgeid eid) {
        //     scoreData[eid] = baseData[eid] + lowerCorrection[eid];
        //     boundData[eid] = baseData[eid] + upperCorrection[eid];
        // });
    }
    ++iterationReached;
}

bool HyperKatzCentrality::checkGlobalConvergence() {
    bool globalConverged = true;
    for (index i = 0; i < levelMatrices.size(); ++i) {
        if (activeLevel[i])
            globalConverged = globalConverged && checkConvergence(i);
    }

    return globalConverged;
}

bool HyperKatzCentrality::checkConvergence(index i) {
    if (activeRanking[i].size() > k) {
        std::partial_sort(
            activeRanking[i].begin(), activeRanking[i].begin() + k, activeRanking[i].end(),
            [&](edgeid eid1, edgeid eid2) { return lowerBound[i][eid1] > lowerBound[i][eid2]; });

        const edgeid kth = activeRanking[i][k - 1];
        activeRanking[i].erase(
            std::remove_if(activeRanking[i].begin() + k, activeRanking[i].end(),
                           [&](edgeid eid) { return areSufficientlyRanked(i, kth, eid); }),
            activeRanking[i].end());
    }

    if (activeRanking[i].size() > k)
        return false;
    if (groupOnly)
        activeLevel[i] = true;
    return true;

    std::sort(activeRanking[i].begin(), activeRanking[i].end(),
              [&](edgeid eid1, edgeid eid2) { return lowerBound[i][eid1] > lowerBound[i][eid2]; });
    for (index j = 1; j < activeRanking[i].size(); ++j) {
        if (!areSufficientlyRanked(i, activeRanking[i][j - 1], activeRanking[i][j]))
            return false;
    }

    activeLevel[i] = true;
    return true;
}

} // namespace NetworKit
