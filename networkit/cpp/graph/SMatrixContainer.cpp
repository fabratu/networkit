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

struct SparsePattern {
    std::vector<index> rowOffsets;
    std::vector<index> columnIndices;
};

} // namespace

template <typename Matrix>
void SMatrixContainer<Matrix>::build(SMatrixType matrixType, VSRMatrix *vsrMatrix,
                                     std::vector<double> *alphaVector, bool patternOnly) {
    reset();

    if (patternOnly) {
        if constexpr (!std::is_same_v<Matrix, ACSRMatrix> && !std::is_same_v<Matrix, DCSRMatrix>)
            throw std::invalid_argument("Pattern-only storage requires ACSRMatrix or DCSRMatrix");
        if (matrixType != SMatrixType::Delta && !std::is_same_v<Matrix, DCSRMatrix>)
            throw std::invalid_argument("Pattern-only level storage requires DCSRMatrix");
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
        const bool parallelPatternBuild =
            patternOnly && matrixType == SMatrixType::Delta
            && (std::is_same_v<Matrix, ACSRMatrix> || std::is_same_v<Matrix, DCSRMatrix>);
        const count numberOfChunks =
            pairs.empty()
                ? 1
                : std::min<count>(Aux::getMaxNumberOfThreads(), static_cast<count>(pairs.size()));
        std::vector<index> chunkBounds(numberOfChunks + 1, 0);
        for (count chunk = 1; chunk < numberOfChunks; ++chunk) {
            index boundary = pairs.size() * chunk / numberOfChunks;
            while (boundary < pairs.size() && pairs[boundary] == pairs[boundary - 1])
                ++boundary;
            chunkBounds[chunk] = boundary;
        }
        chunkBounds[numberOfChunks] = pairs.size();

        std::vector<count> bucketSizes;
        std::vector<std::vector<count>> chunkBucketSizes(numberOfChunks);
        count chunk = 0;
        for (index begin = 0; begin < pairs.size();) {
            index end = begin + 1;
            while (end < pairs.size() && pairs[end] == pairs[begin])
                ++end;
            const count intersectionSize = end - begin;
            maximumIntersection = std::max(maximumIntersection, intersectionSize);
            if (bucketSizes.size() <= intersectionSize)
                bucketSizes.resize(intersectionSize + 1, 0);
            ++bucketSizes[intersectionSize];
            if (parallelPatternBuild) {
                while (chunk + 1 < numberOfChunks && begin >= chunkBounds[chunk + 1])
                    ++chunk;
                if (chunkBucketSizes[chunk].size() <= intersectionSize)
                    chunkBucketSizes[chunk].resize(intersectionSize + 1, 0);
                ++chunkBucketSizes[chunk][intersectionSize];
            }
            begin = end;
        }

        deltaEntries.resize(maximumIntersection + 2);
        if (parallelPatternBuild) {
            std::vector<std::vector<index>> writeOffsets(
                numberOfChunks, std::vector<index>(maximumIntersection + 1, 0));
            for (count s = 1; s <= maximumIntersection; ++s) {
                deltaEntries[s].resize(2 * bucketSizes[s]);
                index offset = 0;
                for (count chunkIndex = 0; chunkIndex < numberOfChunks; ++chunkIndex) {
                    writeOffsets[chunkIndex][s] = offset;
                    if (chunkBucketSizes[chunkIndex].size() > s)
                        offset += 2 * chunkBucketSizes[chunkIndex][s];
                }
            }

#pragma omp parallel for schedule(static)
            for (omp_index chunkIndex = 0; chunkIndex < static_cast<omp_index>(numberOfChunks);
                 ++chunkIndex) {
                auto positions = writeOffsets[chunkIndex];
                for (index begin = chunkBounds[chunkIndex]; begin < chunkBounds[chunkIndex + 1];) {
                    index end = begin + 1;
                    while (end < chunkBounds[chunkIndex + 1] && pairs[end] == pairs[begin])
                        ++end;

                    const count intersectionSize = end - begin;
                    const edgeid eid1 = firstOf(pairs[begin]);
                    const edgeid eid2 = secondOf(pairs[begin]);
                    index &position = positions[intersectionSize];
                    deltaEntries[intersectionSize][position++] = {eid1, eid2};
                    deltaEntries[intersectionSize][position++] = {eid2, eid1};
                    begin = end;
                }
            }
            return;
        }

        for (count s = 1; s <= maximumIntersection; ++s)
            deltaEntries[s].reserve(2 * bucketSizes[s]);
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

    const bool deriveMetadataFromPatterns =
        patternOnly && matrixType == SMatrixType::Delta
        && (std::is_same_v<Matrix, ACSRMatrix> || std::is_same_v<Matrix, DCSRMatrix>);
    if (alphaVector != nullptr && !deriveMetadataFromPatterns) {
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
        return SparsePattern{std::move(rowOffsets), std::move(columnIndices)};
    };

    const auto makeCsrPatternParallel = [&](const std::vector<SparseEntry> &entries) {
        // Four partitions keep the per-partition row-position table cache-friendly. On large
        // matrices, using every worker increases this table by hundreds of MiB and is slower.
        constexpr count maximumGroupingChunks = 4;
        const count numberOfChunks = std::min<count>(
            {maximumGroupingChunks, static_cast<count>(Aux::getMaxNumberOfThreads()),
             static_cast<count>(entries.size())});
        const bool sizeFits =
            dimension == 0 || numberOfChunks <= std::numeric_limits<count>::max() / dimension;
        const bool denseEnough = dimension != 0 && entries.size() / dimension >= 4;
        if (numberOfChunks < 2 || !sizeFits || !denseEnough)
            return makeCsrPattern(entries);

        const auto groupEntries = [&](auto positionType) {
            using Position = decltype(positionType);
            std::vector<Position> chunkRowPositions(numberOfChunks * dimension, 0);
#pragma omp parallel for schedule(static)
            for (omp_index chunk = 0; chunk < static_cast<omp_index>(numberOfChunks); ++chunk) {
                Position *const rowCounts = chunkRowPositions.data() + chunk * dimension;
                const index begin = entries.size() * chunk / numberOfChunks;
                const index end = entries.size() * (chunk + 1) / numberOfChunks;
                for (index entry = begin; entry < end; ++entry)
                    ++rowCounts[entries[entry].row];
            }

            std::vector<index> rowOffsets(dimension + 1, 0);
#pragma omp parallel for schedule(static)
            for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                count rowCount = 0;
                for (count chunk = 0; chunk < numberOfChunks; ++chunk)
                    rowCount += chunkRowPositions[chunk * dimension + row];
                rowOffsets[row + 1] = rowCount;
            }
            std::partial_sum(rowOffsets.begin(), rowOffsets.end(), rowOffsets.begin());

#pragma omp parallel for schedule(static)
            for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                Position position = rowOffsets[row];
                for (count chunk = 0; chunk < numberOfChunks; ++chunk) {
                    Position &chunkPosition = chunkRowPositions[chunk * dimension + row];
                    const Position chunkCount = chunkPosition;
                    chunkPosition = position;
                    position += chunkCount;
                }
            }

            std::vector<index> columnIndices(entries.size());
#pragma omp parallel for schedule(static)
            for (omp_index chunk = 0; chunk < static_cast<omp_index>(numberOfChunks); ++chunk) {
                Position *const positions = chunkRowPositions.data() + chunk * dimension;
                const index begin = entries.size() * chunk / numberOfChunks;
                const index end = entries.size() * (chunk + 1) / numberOfChunks;
                for (index entry = begin; entry < end; ++entry)
                    columnIndices[positions[entries[entry].row]++] = entries[entry].column;
            }

            return SparsePattern{std::move(rowOffsets), std::move(columnIndices)};
        };

        // Every chunk position is an offset into entries, so a 32-bit table is sufficient
        // whenever the complete pattern fits in 32 bits.
        if (entries.size() <= std::numeric_limits<std::uint32_t>::max())
            return groupEntries(std::uint32_t{});
        return groupEntries(index{});
    };

    levelToMatrix.reserve(maximumIntersection + 1);
    if (matrixType == SMatrixType::Level) {
        if constexpr (std::is_same_v<Matrix, ACSRMatrix>) {
            throw std::invalid_argument("ACSRMatrix storage is only supported for delta matrices");
        } else {
            matrices.reserve(maximumIntersection + 1);
            levelToMatrix.resize(maximumIntersection + 1);

            // Group every exact-intersection bucket by row once. Cumulative level matrices can
            // then be populated by independent row ranges instead of scattering entries through
            // shared row positions.
            std::vector<SparsePattern> exactPatterns(maximumIntersection + 1);
            std::vector<bool> hasExactEntries(maximumIntersection + 1, false);
            for (count s = 1; s <= maximumIntersection; ++s) {
                hasExactEntries[s] = !deltaEntries[s].empty();
                exactPatterns[s] = makeCsrPattern(deltaEntries[s]);
                std::vector<SparseEntry>{}.swap(deltaEntries[s]);
            }

            // The level beyond the maximum intersection is always empty. Build the remaining
            // levels from low to high so each exact pattern can be released after its final use.
            // Empty exact levels reuse the preceding cumulative matrix.
            matrices.emplace_back(dimension);
            constexpr index emptyMatrixIndex = 0;
            levelToMatrix[maximumIntersection] = emptyMatrixIndex;

            std::vector<count> rowCounts(dimension, 0);
            count numberOfEntries = 0;
            for (count s = 1; s <= maximumIntersection; ++s)
                numberOfEntries += exactPatterns[s].columnIndices.size();

#pragma omp parallel for schedule(static)
            for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                count entriesInRow = 0;
                for (count s = 1; s <= maximumIntersection; ++s) {
                    entriesInRow +=
                        exactPatterns[s].rowOffsets[row + 1] - exactPatterns[s].rowOffsets[row];
                }
                rowCounts[row] = entriesInRow;
            }

            index currentMatrix = emptyMatrixIndex;
            for (count s = 1; s <= maximumIntersection; ++s) {
                if (s == 1 || hasExactEntries[s - 1]) {
                    std::vector<index> rowOffsets(dimension + 1, 0);
                    for (index row = 0; row < dimension; ++row)
                        rowOffsets[row + 1] = rowOffsets[row] + rowCounts[row];

                    std::vector<index> columnIndices(numberOfEntries);
                    std::vector<double> values;
                    if (!patternOnly)
                        values.resize(numberOfEntries);
#pragma omp parallel for schedule(guided)
                    for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                        index destination = rowOffsets[row];
                        for (count level = s; level <= maximumIntersection; ++level) {
                            const auto &pattern = exactPatterns[level];
                            for (index entry = pattern.rowOffsets[row];
                                 entry < pattern.rowOffsets[row + 1]; ++entry) {
                                columnIndices[destination] = pattern.columnIndices[entry];
                                if (!patternOnly)
                                    values[destination] = 1.0;
                                ++destination;
                            }
                        }
                    }

                    if constexpr (std::is_same_v<Matrix, CSRMatrix>) {
                        matrices.emplace_back(dimension, dimension, std::move(rowOffsets),
                                              std::move(columnIndices), std::move(values), 0.0,
                                              false);
                    } else {
                        if (patternOnly) {
                            matrices.emplace_back(dimension, dimension, std::move(rowOffsets),
                                                  std::move(columnIndices), false);
                        } else {
                            matrices.emplace_back(dimension, dimension, rowOffsets, columnIndices,
                                                  values, 0.0, false);
                        }
                    }
                    currentMatrix = matrices.size() - 1;
                }
                levelToMatrix[s - 1] = currentMatrix;

#pragma omp parallel for schedule(static)
                for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                    rowCounts[row] -=
                        exactPatterns[s].rowOffsets[row + 1] - exactPatterns[s].rowOffsets[row];
                }
                numberOfEntries -= exactPatterns[s].columnIndices.size();
                exactPatterns[s] = {};
            }
        }
    } else {
        matrices.reserve(maximumIntersection + 1);

        if constexpr (std::is_same_v<Matrix, ACSRMatrix>) {
            if (patternOnly) {
                combinedPatternBuilt = true;
                std::vector<SparsePattern> exactPatterns(maximumIntersection + 2);
                count numberOfEntries = 0;
                for (count s = 1; s <= maximumIntersection; ++s) {
                    exactPatterns[s] = makeCsrPatternParallel(deltaEntries[s]);
                    numberOfEntries += exactPatterns[s].columnIndices.size();
                    std::vector<SparseEntry>{}.swap(deltaEntries[s]);
                }
                exactPatterns[maximumIntersection + 1].rowOffsets.assign(dimension + 1, 0);

                if (alphaVector != nullptr)
                    alphaVector->assign(maximumIntersection + 1, 1.0);
                if (alphaVector != nullptr || !maximumLevels.empty()) {
                    std::vector<count> degrees(dimension, 0);
                    for (count s = maximumIntersection; s > 0; --s) {
                        count maximumDegree = 0;
                        const auto &pattern = exactPatterns[s];
#pragma omp parallel for reduction(max : maximumDegree) schedule(static)
                        for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                            const count exactDegree =
                                pattern.rowOffsets[row + 1] - pattern.rowOffsets[row];
                            degrees[row] += exactDegree;
                            maximumDegree = std::max(maximumDegree, degrees[row]);
                            if (!maximumLevels.empty() && exactDegree != 0)
                                maximumLevels[row] = std::max(maximumLevels[row], s);
                        }
                        if (alphaVector != nullptr) {
                            (*alphaVector)[s - 1] =
                                1.0 / (static_cast<double>(maximumDegree) + 1.0);
                        }
                    }
                }

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
#pragma omp parallel for schedule(static)
                    for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                        count rowEntries = 0;
                        for (count s = 1; s <= maximumIntersection; ++s) {
                            rowEntries += exactPatterns[s].rowOffsets[row + 1]
                                          - exactPatterns[s].rowOffsets[row];
                        }
                        combinedRowIdx[row + 1] = rowEntries;
                    }
                    std::partial_sum(combinedRowIdx.begin(), combinedRowIdx.end(),
                                     combinedRowIdx.begin());

                    combinedColumnIdx.resize(numberOfEntries);
                    combinedLevels.resize(numberOfEntries);
#pragma omp parallel for schedule(guided)
                    for (omp_index row = 0; row < static_cast<omp_index>(dimension); ++row) {
                        index destination = combinedRowIdx[row];
                        for (count s = 1; s <= maximumIntersection; ++s) {
                            const auto &pattern = exactPatterns[s];
                            for (index entry = pattern.rowOffsets[row];
                                 entry < pattern.rowOffsets[row + 1]; ++entry) {
                                combinedColumnIdx[destination] = pattern.columnIndices[entry];
                                combinedLevels[destination] = s;
                                ++destination;
                            }
                        }
                    }
                }

                for (count s = 1; s <= maximumIntersection + 1; ++s) {
                    auto &pattern = exactPatterns[s];
                    matrices.emplace_back(dimension, dimension, s, std::move(pattern.rowOffsets),
                                          std::move(pattern.columnIndices), false, true);
                    levelToMatrix.push_back(matrices.size() - 1);
                }
            } else {
                for (count s = 1; s <= maximumIntersection + 1; ++s) {
                    auto [rowOffsets, columnIndices] = makeCsrPattern(deltaEntries[s]);
                    std::vector<double> values(columnIndices.size() * s, 1.0);
                    matrices.emplace_back(dimension, dimension, s, std::move(rowOffsets),
                                          std::move(columnIndices), std::move(values), false);
                    std::vector<SparseEntry>{}.swap(deltaEntries[s]);
                    levelToMatrix.push_back(matrices.size() - 1);
                }
            }
        } else {
            if (patternOnly) {
                if constexpr (std::is_same_v<Matrix, DCSRMatrix>) {
                    std::vector<SparsePattern> exactPatterns(maximumIntersection + 2);
                    for (count s = 1; s <= maximumIntersection; ++s) {
                        exactPatterns[s] = makeCsrPatternParallel(deltaEntries[s]);
                        std::vector<SparseEntry>{}.swap(deltaEntries[s]);
                    }
                    exactPatterns[maximumIntersection + 1].rowOffsets.assign(dimension + 1, 0);

                    if (alphaVector != nullptr)
                        alphaVector->assign(maximumIntersection + 1, 1.0);
                    if (alphaVector != nullptr || !maximumLevels.empty()) {
                        std::vector<count> degrees(dimension, 0);
                        for (count s = maximumIntersection; s > 0; --s) {
                            count maximumDegree = 0;
                            const auto &pattern = exactPatterns[s];
#pragma omp parallel for reduction(max : maximumDegree) schedule(static)
                            for (omp_index row = 0; row < static_cast<omp_index>(dimension);
                                 ++row) {
                                const count exactDegree =
                                    pattern.rowOffsets[row + 1] - pattern.rowOffsets[row];
                                degrees[row] += exactDegree;
                                maximumDegree = std::max(maximumDegree, degrees[row]);
                                if (!maximumLevels.empty() && exactDegree != 0)
                                    maximumLevels[row] = std::max(maximumLevels[row], s);
                            }
                            if (alphaVector != nullptr)
                                (*alphaVector)[s - 1] =
                                    1.0 / (static_cast<double>(maximumDegree) + 1.0);
                        }
                    }

                    for (count s = 1; s <= maximumIntersection + 1; ++s) {
                        auto &pattern = exactPatterns[s];
                        matrices.emplace_back(dimension, dimension, std::move(pattern.rowOffsets),
                                              std::move(pattern.columnIndices), false);
                        levelToMatrix.push_back(matrices.size() - 1);
                    }
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
                                              std::move(columnIndices), std::move(values), 0.0,
                                              false);
                    } else {
                        matrices.emplace_back(dimension, dimension, rowOffsets, columnIndices,
                                              values, 0.0, false);
                    }
                    std::vector<SparseEntry>{}.swap(deltaEntries[s]);
                    levelToMatrix.push_back(matrices.size() - 1);
                }
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
                if (matrix.k == 1) {
                    if (matrix.compressedRows) {
#pragma omp for schedule(static)
                        for (omp_index compressedRow = 0;
                             compressedRow < static_cast<omp_index>(matrix.nonZeroRows.size());
                             ++compressedRow) {
                            const index row = matrix.nonZeroRows[compressedRow];
                            double sum = output(row, 0);
                            for (index entry = matrix.rowIdx[compressedRow];
                                 entry < matrix.rowIdx[compressedRow + 1]; ++entry)
                                sum += input(matrix.columnIdx[entry], 0);
                            output(row, 0) = sum;
                        }
                        continue;
                    }
#pragma omp for schedule(static)
                    for (omp_index row = 0; row < static_cast<omp_index>(hGraph.upperEdgeIdBound());
                         ++row) {
                        double sum = output(row, 0);
                        for (index entry = matrix.rowIdx[row]; entry < matrix.rowIdx[row + 1];
                             ++entry)
                            sum += input(matrix.columnIdx[entry], 0);
                        output(row, 0) = sum;
                    }
                    continue;
                }
                if (matrix.compressedRows) {
#pragma omp for schedule(guided)
                    for (omp_index compressedRow = 0;
                         compressedRow < static_cast<omp_index>(matrix.nonZeroRows.size());
                         ++compressedRow) {
                        double *const outputRow = &output(matrix.nonZeroRows[compressedRow], 0);
                        for (index entry = matrix.rowIdx[compressedRow];
                             entry < matrix.rowIdx[compressedRow + 1]; ++entry) {
                            const double *const inputRow = &input(matrix.columnIdx[entry], 0);
                            std::transform(inputRow, inputRow + matrix.k, outputRow, outputRow,
                                           std::plus<double>{});
                        }
                    }
                    continue;
                }
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
