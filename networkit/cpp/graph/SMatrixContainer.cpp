#include <algorithm>
#include <stdexcept>
#include <unordered_map>

#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

void SMatrixContainer::build(SMatrixType matrixType) {
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

    levelToMatrix.reserve(maximumIntersection + 1);
    if (matrixType == SMatrixType::Level) {
        // A level matrix changes from s - 1 to s exactly when a pair has intersection size s - 1.
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

                std::sort(
                    triplets.begin(), triplets.end(), [](const Triplet &lhs, const Triplet &rhs) {
                        return lhs.row < rhs.row || (lhs.row == rhs.row && lhs.column < rhs.column);
                    });
                matrices.emplace_back(dimension, triplets, 0.0, true);
            }

            levelToMatrix.push_back(matrices.size() - 1);
        }
    } else {
        // All delta levels without an exact interaction share this matrix.
        matrices.emplace_back(dimension);
        constexpr index emptyMatrixIndex = 0;

        for (count s = 1; s <= maximumIntersection + 1; ++s) {
            std::vector<Triplet> triplets;
            for (edgeid eid1 = 0; eid1 < intersectionSizes.size(); ++eid1) {
                for (const auto &[eid2, intersectionSize] : intersectionSizes[eid1]) {
                    if (intersectionSize == s) {
                        triplets.push_back({eid1, eid2, 1.0});
                        triplets.push_back({eid2, eid1, 1.0});
                    }
                }
            }

            if (triplets.empty()) {
                levelToMatrix.push_back(emptyMatrixIndex);
                continue;
            }

            std::sort(triplets.begin(), triplets.end(), [](const Triplet &lhs, const Triplet &rhs) {
                return lhs.row < rhs.row || (lhs.row == rhs.row && lhs.column < rhs.column);
            });
            matrices.emplace_back(dimension, triplets, 0.0, true);
            levelToMatrix.push_back(matrices.size() - 1);
        }
    }

    type = matrixType;
    built = true;
}

void SMatrixContainer::reset() noexcept {
    matrices.clear();
    levelToMatrix.clear();
    built = false;
}

const CSRMatrix &SMatrixContainer::getMatrix(count s) const {
    if (s == 0)
        throw std::invalid_argument("The s-level must be positive");

    if (!built)
        throw std::runtime_error("SMatrixContainer has not been built");

    if (s > getMaxLevel())
        throw std::out_of_range("The s-level exceeds the maximum level");

    return matrices[levelToMatrix[s - 1]];
}

SMatrixType SMatrixContainer::getType() const {
    if (!built)
        throw std::runtime_error("SMatrixContainer has not been built");

    return type;
}

} // namespace NetworKit
