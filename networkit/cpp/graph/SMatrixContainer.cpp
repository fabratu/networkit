#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

template <typename Matrix>
void SMatrixContainer<Matrix>::build(SMatrixType matrixType, VSRMatrix *vsrMatrix,
                                     std::vector<double> *alphaVector) {
    reset();

    const count dimension = hGraph.upperEdgeIdBound();
    std::vector<std::unordered_map<edgeid, count>> intersectionSizes(dimension);
    count maximumIntersection = 0;

    // Count every hyperedge intersection in one pass over the node incidences. Pairs are stored
    // only for the smaller edge id, so each common node contributes exactly once.
    hGraph.forNodes([&](node u) {
        const auto &incidentEdges = hGraph.edgesOf(u);
        for (auto first = incidentEdges.begin(); first != incidentEdges.end(); ++first) {
            for (auto second = std::next(first); second != incidentEdges.end(); ++second) {
                const auto [eid1, eid2] = std::minmax(*first, *second);
                auto &intersectionSize = intersectionSizes[eid1][eid2];
                maximumIntersection = std::max(maximumIntersection, ++intersectionSize);
            }
        }
    });

    // Delta entries are bucketed by their exact intersection level. Besides avoiding a complete
    // scan of all intersections for every level, these directed entries also allow all maximum
    // level degrees to be computed without constructing level matrices.
    std::vector<std::vector<Triplet>> deltaEntries;
    if (matrixType == SMatrixType::Delta || alphaVector != nullptr)
        deltaEntries.resize(maximumIntersection + 2);

    std::vector<count> maximumLevels;
    if (vsrMatrix != nullptr)
        maximumLevels.assign(dimension, 1);

    if (!deltaEntries.empty() || !maximumLevels.empty()) {
        for (edgeid eid1 = 0; eid1 < intersectionSizes.size(); ++eid1) {
            for (const auto &[eid2, intersectionSize] : intersectionSizes[eid1]) {
                if (!deltaEntries.empty()) {
                    deltaEntries[intersectionSize].push_back({eid1, eid2, 1.0});
                    deltaEntries[intersectionSize].push_back({eid2, eid1, 1.0});
                }
                if (!maximumLevels.empty()) {
                    maximumLevels[eid1] = std::max(maximumLevels[eid1], intersectionSize);
                    maximumLevels[eid2] = std::max(maximumLevels[eid2], intersectionSize);
                }
            }
        }
    }

    if (alphaVector != nullptr) {
        alphaVector->assign(maximumIntersection + 1, 1.0);
        std::vector<count> degrees(dimension, 0);
        count maximumDegree = 0;
        for (count s = maximumIntersection; s > 0; --s) {
            for (const auto &entry : deltaEntries[s])
                maximumDegree = std::max(maximumDegree, ++degrees[entry.row]);
            (*alphaVector)[s - 1] = 1.0 / (static_cast<double>(maximumDegree) + 1.0);
        }
    }

    levelToMatrix.reserve(maximumIntersection + 1);
    if (matrixType == SMatrixType::Level) {
        if constexpr (std::is_same_v<Matrix, ACSRMatrix>) {
            throw std::invalid_argument("ACSRMatrix storage is only supported for delta matrices");
        } else {
            // A level matrix changes from s - 1 to s exactly when a pair has intersection size
            // s - 1.
            std::vector<bool> changesAtLevel(maximumIntersection + 2, false);
            changesAtLevel[1] = true;
            changesAtLevel[maximumIntersection + 1] = true;
            for (const auto &row : intersectionSizes) {
                for (const auto &entry : row)
                    changesAtLevel[entry.second + 1] = true;
            }

            for (count s = 1; s <= maximumIntersection + 1; ++s) {
                if (changesAtLevel[s]) {
                    std::vector<Triplet> triplets;
                    for (edgeid eid1 = 0; eid1 < intersectionSizes.size(); ++eid1) {
                        for (const auto &[eid2, intersectionSize] : intersectionSizes[eid1]) {
                            if (intersectionSize >= s) {
                                triplets.push_back({eid1, eid2, 1.0});
                                triplets.push_back({eid2, eid1, 1.0});
                            }
                        }
                    }

                    std::sort(triplets.begin(), triplets.end(),
                              [](const Triplet &lhs, const Triplet &rhs) {
                                  return lhs.row < rhs.row
                                         || (lhs.row == rhs.row && lhs.column < rhs.column);
                              });
                    matrices.emplace_back(dimension, triplets, 0.0, true);
                }

                levelToMatrix.push_back(matrices.size() - 1);
            }
        }
    } else {
        // All information needed by delta matrices has been extracted into compact level buckets.
        // Release the much heavier hash-map representation before allocating matrix values.
        intersectionSizes.clear();
        intersectionSizes.shrink_to_fit();
        matrices.reserve(maximumIntersection + 1);

        if constexpr (std::is_same_v<Matrix, ACSRMatrix>) {
            for (count s = 1; s <= maximumIntersection + 1; ++s) {
                matrices.emplace_back(dimension, dimension, s, deltaEntries[s], 1.0, false);
                std::vector<Triplet>{}.swap(deltaEntries[s]);
                levelToMatrix.push_back(matrices.size() - 1);
            }
        } else {
            // All delta levels without an exact interaction share this matrix.
            matrices.emplace_back(dimension);
            constexpr index emptyMatrixIndex = 0;

            for (count s = 1; s <= maximumIntersection + 1; ++s) {
                if (deltaEntries[s].empty()) {
                    levelToMatrix.push_back(emptyMatrixIndex);
                    continue;
                }

                matrices.emplace_back(dimension, deltaEntries[s], 0.0, false);
                std::vector<Triplet>{}.swap(deltaEntries[s]);
                levelToMatrix.push_back(matrices.size() - 1);
            }
        }
    }

    if (vsrMatrix != nullptr) {
        std::vector<count> rowLengths;
        rowLengths.reserve(hGraph.numberOfEdges());
        hGraph.forEdges([&](edgeid eid) { rowLengths.push_back(maximumLevels[eid]); });

        const count maximumRowLength =
            rowLengths.empty() ? 1 : *std::max_element(rowLengths.begin(), rowLengths.end());
        *vsrMatrix =
            VSRMatrix(hGraph.numberOfEdges(), maximumRowLength + 1, maximumRowLength, rowLengths);
    }

    type = matrixType;
    built = true;
}

template <typename Matrix>
void SMatrixContainer<Matrix>::reset() noexcept {
    matrices.clear();
    levelToMatrix.clear();
    built = false;
}

template <typename Matrix>
const Matrix &SMatrixContainer<Matrix>::getMatrix(count s) const {
    if (s == 0)
        throw std::invalid_argument("The s-level must be positive");

    if (!built)
        throw std::runtime_error("SMatrixContainer has not been built");

    if (s > getMaxLevel())
        throw std::out_of_range("The s-level exceeds the maximum level");

    return matrices[levelToMatrix[s - 1]];
}

template <typename Matrix>
SMatrixType SMatrixContainer<Matrix>::getType() const {
    if (!built)
        throw std::runtime_error("SMatrixContainer has not been built");

    return type;
}

template class SMatrixContainer<CSRMatrix>;
template class SMatrixContainer<DCSRMatrix>;
template class SMatrixContainer<ACSRMatrix>;

} // namespace NetworKit
