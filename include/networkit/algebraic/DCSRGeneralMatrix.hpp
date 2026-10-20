/*
 * DCSRGeneralMatrix.hpp
 *
 * A doubly compressed sparse row matrix implementation.
 */

#ifndef NETWORKIT_ALGEBRAIC_DCSR_GENERAL_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_DCSR_GENERAL_MATRIX_HPP_

#include <algorithm>
#include <cassert>
#include <numeric>
#include <utility>
#include <vector>

#include <omp.h>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Sparse matrix stored in doubly compressed sparse row (DCSR) format.
 *
 * Unlike CSR, DCSR stores row offsets only for non-empty rows. This makes it
 * suitable for matrices with a large number of empty rows. The supported
 * operations are intentionally limited to matrix-vector multiplication and
 * basic matrix properties.
 */
template <class ValueType>
class DCSRGeneralMatrix {
    std::vector<index> nonZeroRows;
    std::vector<index> rowIdx;
    std::vector<index> columnIdx;
    std::vector<ValueType> nonZeros;

    count nRows;
    count nCols;
    bool isSorted;
    ValueType zero;
    bool unitValues = false;

public:
    /** Default constructor. */
    DCSRGeneralMatrix()
        : nonZeroRows(), rowIdx(1, 0), columnIdx(), nonZeros(), nRows(0), nCols(0), isSorted(true),
          zero(0) {}

    /** Constructs an empty square matrix. */
    DCSRGeneralMatrix(count dimension, ValueType zero = 0)
        : nonZeroRows(), rowIdx(1, 0), columnIdx(), nonZeros(), nRows(dimension), nCols(dimension),
          isSorted(true), zero(zero) {}

    /** Constructs an empty rectangular matrix. */
    DCSRGeneralMatrix(count nRows, count nCols, ValueType zero = 0)
        : nonZeroRows(), rowIdx(1, 0), columnIdx(), nonZeros(), nRows(nRows), nCols(nCols),
          isSorted(true), zero(zero) {}

    /** Constructs a square matrix from triplets. */
    DCSRGeneralMatrix(count dimension, const std::vector<Triplet> &triplets, ValueType zero = 0,
                      bool isSorted = false)
        : DCSRGeneralMatrix(dimension, dimension, triplets, zero, isSorted) {}

    /** Constructs a rectangular matrix from triplets. */
    DCSRGeneralMatrix(count nRows, count nCols, const std::vector<Triplet> &triplets,
                      ValueType zero = 0, bool isSorted = false)
        : nonZeroRows(), rowIdx(), columnIdx(triplets.size()), nonZeros(triplets.size()),
          nRows(nRows), nCols(nCols), isSorted(isSorted), zero(zero) {
        std::vector<index> order(triplets.size());
        std::iota(order.begin(), order.end(), index{0});
        for (const auto &triplet : triplets) {
            assert(triplet.row < nRows);
            assert(triplet.column < nCols);
        }
        std::stable_sort(order.begin(), order.end(), [&triplets](index lhs, index rhs) {
            return triplets[lhs].row < triplets[rhs].row;
        });

        nonZeroRows.reserve(std::min<count>(nRows, triplets.size()));
        rowIdx.reserve(nonZeroRows.capacity() + 1);
        rowIdx.push_back(0);
        for (index i = 0; i < order.size(); ++i) {
            const auto &triplet = triplets[order[i]];
            if (i == 0 || triplet.row != nonZeroRows.back()) {
                if (i > 0)
                    rowIdx.push_back(i);
                nonZeroRows.push_back(triplet.row);
            }
            columnIdx[i] = triplet.column;
            nonZeros[i] = triplet.value;
        }
        if (!triplets.empty())
            rowIdx.push_back(triplets.size());
    }

    /** Constructs a matrix from columns and values grouped by row. */
    DCSRGeneralMatrix(count nRows, count nCols, const std::vector<std::vector<index>> &columnIdx,
                      const std::vector<std::vector<ValueType>> &values, ValueType zero = 0,
                      bool isSorted = false)
        : nonZeroRows(), rowIdx(1, 0), columnIdx(), nonZeros(), nRows(nRows), nCols(nCols),
          isSorted(isSorted), zero(zero) {
        assert(columnIdx.size() == nRows);
        assert(values.size() == nRows);

        for (index i = 0; i < nRows; ++i) {
            assert(columnIdx[i].size() == values[i].size());
            if (columnIdx[i].empty())
                continue;

            nonZeroRows.push_back(i);
            this->columnIdx.insert(this->columnIdx.end(), columnIdx[i].begin(), columnIdx[i].end());
            nonZeros.insert(nonZeros.end(), values[i].begin(), values[i].end());
            rowIdx.push_back(nonZeros.size());
        }
    }

    /** Constructs a matrix from CSR arrays, compressing its empty rows. */
    DCSRGeneralMatrix(count nRows, count nCols, const std::vector<index> &rowIdx,
                      const std::vector<index> &columnIdx, const std::vector<ValueType> &nonZeros,
                      ValueType zero = 0, bool isSorted = false)
        : nonZeroRows(), rowIdx(1, 0), columnIdx(columnIdx), nonZeros(nonZeros), nRows(nRows),
          nCols(nCols), isSorted(isSorted), zero(zero) {
        assert(rowIdx.size() == nRows + 1);
        assert(columnIdx.size() == nonZeros.size());
        assert(rowIdx.back() == nonZeros.size());

        nonZeroRows.reserve(std::min<count>(nRows, nonZeros.size()));
        this->rowIdx.reserve(nonZeroRows.capacity() + 1);
        for (index i = 0; i < nRows; ++i) {
            if (rowIdx[i] != rowIdx[i + 1]) {
                nonZeroRows.push_back(i);
                this->rowIdx.push_back(rowIdx[i + 1]);
            }
        }
    }

    /** Takes ownership of a CSR pattern and stores its entries as implicit unit values. */
    DCSRGeneralMatrix(count nRows, count nCols, std::vector<index> &&rowOffsets,
                      std::vector<index> &&columnIndices, bool isSorted = false)
        : nonZeroRows(), rowIdx(1, 0), columnIdx(std::move(columnIndices)), nonZeros(),
          nRows(nRows), nCols(nCols), isSorted(isSorted), zero(0), unitValues(true) {
        assert(rowOffsets.size() == nRows + 1);
        assert(rowOffsets.front() == 0 && rowOffsets.back() == columnIdx.size());
        nonZeroRows.reserve(std::min<count>(nRows, columnIdx.size()));
        this->rowIdx.reserve(nonZeroRows.capacity() + 1);
        for (index row = 0; row < nRows; ++row) {
            if (rowOffsets[row] != rowOffsets[row + 1]) {
                nonZeroRows.push_back(row);
                this->rowIdx.push_back(rowOffsets[row + 1]);
            }
        }
    }

    DCSRGeneralMatrix(const DCSRGeneralMatrix &other) = default;
    DCSRGeneralMatrix(DCSRGeneralMatrix &&other) noexcept = default;
    ~DCSRGeneralMatrix() = default;

    DCSRGeneralMatrix &operator=(const DCSRGeneralMatrix &other) = default;
    DCSRGeneralMatrix &operator=(DCSRGeneralMatrix &&other) noexcept = default;

    /** @return Number of rows. */
    count numberOfRows() const noexcept { return nRows; }

    /** @return Number of columns. */
    count numberOfColumns() const noexcept { return nCols; }

    /** @return The matrix zero element. */
    ValueType getZero() const noexcept { return zero; }

    /** @return Number of stored entries. */
    count nnz() const noexcept { return columnIdx.size(); }

    /** Whether stored entries have implicit unit values. */
    bool hasUnitValues() const noexcept { return unitValues; }

    /** @return Whether column indices within each row are sorted. */
    bool sorted() const noexcept { return isSorted; }

    /** @return Number of stored entries in row @a i. */
    count nnzInRow(index i) const {
        assert(i < nRows);
        const auto it = std::lower_bound(nonZeroRows.begin(), nonZeroRows.end(), i);
        if (it == nonZeroRows.end() || *it != i)
            return 0;
        const index compressedRow = it - nonZeroRows.begin();
        return rowIdx[compressedRow + 1] - rowIdx[compressedRow];
    }

    /** @return Largest number of stored entries in any row. */
    count maximumRowDegree() const noexcept {
        count maximum = 0;
        for (index row = 0; row < nonZeroRows.size(); ++row)
            maximum = std::max(maximum, rowIdx[row + 1] - rowIdx[row]);
        return maximum;
    }

    /** Multiplies this matrix with a column vector. */
    Vector operator*(const Vector &vector) const {
        assert(!vector.isTransposed());
        assert(nCols == vector.getDimension());

        Vector result(nRows, zero);
        if (unitValues) {
#pragma omp parallel for
            for (omp_index r = 0; r < static_cast<omp_index>(nonZeroRows.size()); ++r) {
                ValueType sum = zero;
                for (index k = rowIdx[r]; k < rowIdx[r + 1]; ++k)
                    sum += vector[columnIdx[k]];
                result[nonZeroRows[r]] = sum;
            }
            return result;
        }
#pragma omp parallel for
        for (omp_index r = 0; r < static_cast<omp_index>(nonZeroRows.size()); ++r) {
            ValueType sum = zero;
            for (index k = rowIdx[r]; k < rowIdx[r + 1]; ++k)
                sum += nonZeros[k] * vector[columnIdx[k]];
            result[nonZeroRows[r]] = sum;
        }

        return result;
    }

    /** Adds the product of an implicit-unit pattern and @a input to @a output. */
    void addUnitValuesInto(const Vector &input, Vector &output) const {
        assert(unitValues);
        assert(!input.isTransposed() && !output.isTransposed());
        assert(input.getDimension() == nCols && output.getDimension() == nRows);
#pragma omp parallel for
        for (omp_index r = 0; r < static_cast<omp_index>(nonZeroRows.size()); ++r) {
            ValueType sum = zero;
            for (index k = rowIdx[r]; k < rowIdx[r + 1]; ++k)
                sum += input[columnIdx[k]];
            output[nonZeroRows[r]] += sum;
        }
    }

    /** Overwrites @a output with the product of the unit-valued sparsity pattern and @a input. */
    void multiplyUnitValuesInto(const Vector &input, Vector &output) const {
        assert(!input.isTransposed() && !output.isTransposed());
        assert(input.getDimension() == nCols && output.getDimension() == nRows);
        output.fill(zero);
#pragma omp parallel for
        for (omp_index r = 0; r < static_cast<omp_index>(nonZeroRows.size()); ++r) {
            ValueType sum = zero;
            for (index k = rowIdx[r]; k < rowIdx[r + 1]; ++k)
                sum += input[columnIdx[k]];
            output[nonZeroRows[r]] = sum;
        }
    }

};

} // namespace NetworKit

#endif // NETWORKIT_ALGEBRAIC_DCSR_GENERAL_MATRIX_HPP_
