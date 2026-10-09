#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <tlx/sort/parallel_mergesort.hpp>

#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/graph/SMatrixContainer.hpp>

namespace NetworKit {

namespace {

struct EdgePair {
    edgeid first;
    edgeid second;

    bool operator==(const EdgePair &) const = default;
};

struct SparseEntry {
    edgeid row;
    edgeid column;
};

} // namespace

template <typename Matrix>
void SMatrixContainer<Matrix>::build(SMatrixType matrixType, VSRMatrix *vsrMatrix,
                                     std::vector<double> *alphaVector, bool patternOnly) {
    reset();

    if (patternOnly) {
        if constexpr (!std::is_same_v<Matrix, ACSRMatrix>)
            throw std::invalid_argument("Pattern-only storage is supported only for ACSRMatrix");
        if (matrixType != SMatrixType::Delta)
            throw std::invalid_argument("Pattern-only ACSRMatrix storage requires delta matrices");
    }

    const count dimension = hGraph.upperEdgeIdBound();
    const count nodeBound = hGraph.upperNodeIdBound();
    std::vector<index> pairOffsets(nodeBound + 1, 0);
    for (node u = 0; u < nodeBound; ++u) {
        if (!hGraph.hasNode(u))
            continue;
        const count degree = hGraph.edgesOf(u).size();
        pairOffsets[u + 1] = degree > 1 ? degree * (degree - 1) / 2 : 0;
    }
    std::partial_sum(pairOffsets.begin(), pairOffsets.end(), pairOffsets.begin());

    count maximumIntersection = 0;
    std::vector<std::vector<SparseEntry>> deltaEntries;
    std::vector<count> maximumLevels;
    if (vsrMatrix != nullptr)
        maximumLevels.assign(dimension, 1);

    const auto createBuckets = [&](const auto &pairs, auto firstOf, auto secondOf) {
        for (index begin = 0; begin < pairs.size();) {
            index end = begin + 1;
            while (end < pairs.size() && pairs[end] == pairs[begin])
                ++end;
            maximumIntersection = std::max(maximumIntersection, end - begin);
            begin = end;
        }

        deltaEntries.resize(maximumIntersection + 2);
        for (index begin = 0; begin < pairs.size();) {
            index end = begin + 1;
            while (end < pairs.size() && pairs[end] == pairs[begin])
                ++end;

            const count intersectionSize = end - begin;
            const edgeid eid1 = firstOf(pairs[begin]);
            const edgeid eid2 = secondOf(pairs[begin]);
            deltaEntries[intersectionSize].push_back({eid1, eid2});
            deltaEntries[intersectionSize].push_back({eid2, eid1});
            if (!maximumLevels.empty()) {
                maximumLevels[eid1] = std::max(maximumLevels[eid1], intersectionSize);
                maximumLevels[eid2] = std::max(maximumLevels[eid2], intersectionSize);
            }
            begin = end;
        }
    };

    if (dimension <= std::numeric_limits<std::uint32_t>::max()) {
        std::vector<std::uint64_t> pairs(pairOffsets.back());
        hGraph.parallelForNodes([&](node u) {
            index position = pairOffsets[u];
            const auto &incidentEdges = hGraph.edgesOf(u);
            for (auto first = incidentEdges.begin(); first != incidentEdges.end(); ++first) {
                for (auto second = std::next(first); second != incidentEdges.end(); ++second) {
                    const auto [eid1, eid2] = std::minmax(*first, *second);
                    pairs[position++] =
                        (static_cast<std::uint64_t>(eid1) << 32) | static_cast<std::uint32_t>(eid2);
                }
            }
        });
        if (pairs.size() > 1) {
            tlx::parallel_mergesort(pairs.begin(), pairs.end(), std::less<std::uint64_t>{},
                                    std::min<count>(pairs.size(), Aux::getMaxNumberOfThreads()));
        }
        createBuckets(
            pairs, [](std::uint64_t pair) { return static_cast<edgeid>(pair >> 32); },
            [](std::uint64_t pair) {
                return static_cast<edgeid>(static_cast<std::uint32_t>(pair));
            });
    } else {
        std::vector<EdgePair> pairs(pairOffsets.back());
        hGraph.parallelForNodes([&](node u) {
            index position = pairOffsets[u];
            const auto &incidentEdges = hGraph.edgesOf(u);
            for (auto first = incidentEdges.begin(); first != incidentEdges.end(); ++first) {
                for (auto second = std::next(first); second != incidentEdges.end(); ++second) {
                    const auto [eid1, eid2] = std::minmax(*first, *second);
                    pairs[position++] = {eid1, eid2};
                }
            }
        });
        const auto less = [](const EdgePair &lhs, const EdgePair &rhs) {
            return lhs.first < rhs.first || (lhs.first == rhs.first && lhs.second < rhs.second);
        };
        if (pairs.size() > 1) {
            tlx::parallel_mergesort(pairs.begin(), pairs.end(), less,
                                    std::min<count>(pairs.size(), Aux::getMaxNumberOfThreads()));
        }
        createBuckets(
            pairs, [](const EdgePair &pair) { return pair.first; },
            [](const EdgePair &pair) { return pair.second; });
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

    const auto makeCsrPattern = [&](const std::vector<SparseEntry> &entries) {
        std::vector<index> rowOffsets(dimension + 1, 0);
        for (const auto &entry : entries)
            ++rowOffsets[entry.row + 1];
        std::partial_sum(rowOffsets.begin(), rowOffsets.end(), rowOffsets.begin());

        std::vector<index> columnIndices(entries.size());
        auto positions = rowOffsets;
        for (const auto &entry : entries)
            columnIndices[positions[entry.row]++] = entry.column;
        return std::make_pair(std::move(rowOffsets), std::move(columnIndices));
    };

    levelToMatrix.reserve(maximumIntersection + 1);
    if (matrixType == SMatrixType::Level) {
        if constexpr (std::is_same_v<Matrix, ACSRMatrix>) {
            throw std::invalid_argument("ACSRMatrix storage is only supported for delta matrices");
        } else {
            matrices.reserve(maximumIntersection + 1);
            levelToMatrix.resize(maximumIntersection + 1);

            // The level beyond the maximum intersection is always empty. Descending from there,
            // add one exact-level bucket at a time. Empty buckets reuse the matrix for the next
            // higher level.
            matrices.emplace_back(dimension);
            index currentMatrix = 0;
            levelToMatrix[maximumIntersection] = currentMatrix;

            std::vector<count> rowCounts(dimension, 0);
            count numberOfEntries = 0;
            for (count s = maximumIntersection; s > 0; --s) {
                if (deltaEntries[s].empty()) {
                    levelToMatrix[s - 1] = currentMatrix;
                    continue;
                }

                for (const auto &entry : deltaEntries[s])
                    ++rowCounts[entry.row];
                numberOfEntries += deltaEntries[s].size();

                std::vector<index> rowOffsets(dimension + 1, 0);
                for (index row = 0; row < dimension; ++row)
                    rowOffsets[row + 1] = rowOffsets[row] + rowCounts[row];

                std::vector<index> columnIndices(numberOfEntries);
                std::vector<double> values(numberOfEntries, 1.0);
                auto positions = rowOffsets;
                for (count level = s; level <= maximumIntersection; ++level) {
                    for (const auto &entry : deltaEntries[level])
                        columnIndices[positions[entry.row]++] = entry.column;
                }

                if constexpr (std::is_same_v<Matrix, CSRMatrix>) {
                    matrices.emplace_back(dimension, dimension, std::move(rowOffsets),
                                          std::move(columnIndices), std::move(values), 0.0, false);
                } else {
                    matrices.emplace_back(dimension, dimension, rowOffsets, columnIndices, values,
                                          0.0, false);
                }
                currentMatrix = matrices.size() - 1;
                levelToMatrix[s - 1] = currentMatrix;
            }
        }
    } else {
        matrices.reserve(maximumIntersection + 1);

        if constexpr (std::is_same_v<Matrix, ACSRMatrix>) {
            if (patternOnly) {
                combinedPatternBuilt = true;
                count numberOfEntries = 0;
                for (count s = 1; s <= maximumIntersection; ++s)
                    numberOfEntries += deltaEntries[s].size();

                const count numberOfLevels = maximumIntersection + 1;
                const bool thresholdFits =
                    dimension == 0
                    || numberOfLevels <= std::numeric_limits<count>::max() / dimension;
                // A row-major variable-level pattern avoids almost all per-level scheduling for
                // sparse inputs. Once there are more entries than row/level slots, keeping entries
                // grouped by their fixed level provides better locality and vectorization.
                useRowMajorCombinedPattern =
                    dimension == 0 ? numberOfEntries == 0
                                   : thresholdFits && numberOfEntries <= dimension * numberOfLevels;
                if (useRowMajorCombinedPattern) {
                    combinedRowIdx.assign(dimension + 1, 0);
                    for (count s = 1; s <= maximumIntersection; ++s) {
                        for (const auto &entry : deltaEntries[s])
                            ++combinedRowIdx[entry.row + 1];
                    }
                    std::partial_sum(combinedRowIdx.begin(), combinedRowIdx.end(),
                                     combinedRowIdx.begin());

                    combinedColumnIdx.resize(numberOfEntries);
                    combinedLevels.resize(numberOfEntries);
                    auto positions = combinedRowIdx;
                    for (count s = 1; s <= maximumIntersection; ++s) {
                        for (const auto &entry : deltaEntries[s]) {
                            const index destination = positions[entry.row]++;
                            combinedColumnIdx[destination] = entry.column;
                            combinedLevels[destination] = s;
                        }
                    }
                }
            }

            for (count s = 1; s <= maximumIntersection + 1; ++s) {
                auto [rowOffsets, columnIndices] = makeCsrPattern(deltaEntries[s]);
                if (patternOnly) {
                    matrices.emplace_back(dimension, dimension, s, std::move(rowOffsets),
                                          std::move(columnIndices), false);
                } else {
                    std::vector<double> values(columnIndices.size() * s, 1.0);
                    matrices.emplace_back(dimension, dimension, s, std::move(rowOffsets),
                                          std::move(columnIndices), std::move(values), false);
                }
                std::vector<SparseEntry>{}.swap(deltaEntries[s]);
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

                auto [rowOffsets, columnIndices] = makeCsrPattern(deltaEntries[s]);
                std::vector<double> values(columnIndices.size(), 1.0);
                if constexpr (std::is_same_v<Matrix, CSRMatrix>) {
                    matrices.emplace_back(dimension, dimension, std::move(rowOffsets),
                                          std::move(columnIndices), std::move(values), 0.0, false);
                } else {
                    matrices.emplace_back(dimension, dimension, rowOffsets, columnIndices, values,
                                          0.0, false);
                }
                std::vector<SparseEntry>{}.swap(deltaEntries[s]);
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
void SMatrixContainer<Matrix>::multiplyInto(const VSRMatrix &input, VSRMatrix &output) const {
    if constexpr (!std::is_same_v<Matrix, ACSRMatrix>) {
        throw std::logic_error("Combined multiplication is supported only for ACSRMatrix");
    } else {
        if (!built || !combinedPatternBuilt)
            throw std::logic_error("No combined ACSR delta pattern has been built");

        if (useRowMajorCombinedPattern) {
#pragma omp parallel for schedule(guided)
            for (omp_index row = 0; row < static_cast<omp_index>(hGraph.upperEdgeIdBound());
                 ++row) {
                double *const outputRow = &output(row, 0);
                for (index entry = combinedRowIdx[row]; entry < combinedRowIdx[row + 1]; ++entry) {
                    const double *const inputRow = &input(combinedColumnIdx[entry], 0);
                    const count level = combinedLevels[entry];
                    for (index value = 0; value < level; ++value)
                        outputRow[value] += inputRow[value];
                }
            }
            return;
        }

#pragma omp parallel
        {
            for (const auto &matrix : matrices) {
#pragma omp for schedule(guided)
                for (omp_index row = 0; row < static_cast<omp_index>(hGraph.upperEdgeIdBound());
                     ++row) {
                    double *const outputRow = &output(row, 0);
                    for (index entry = matrix.rowIdx[row]; entry < matrix.rowIdx[row + 1];
                         ++entry) {
                        const double *const inputRow = &input(matrix.columnIdx[entry], 0);
                        std::transform(inputRow, inputRow + matrix.k, outputRow, outputRow,
                                       std::plus<double>{});
                    }
                }
            }
        }
    }
}

template <typename Matrix>
void SMatrixContainer<Matrix>::reset() noexcept {
    matrices.clear();
    levelToMatrix.clear();
    combinedRowIdx.clear();
    combinedColumnIdx.clear();
    combinedLevels.clear();
    built = false;
    combinedPatternBuilt = false;
    useRowMajorCombinedPattern = false;
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
