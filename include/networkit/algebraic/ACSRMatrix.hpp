#ifndef NETWORKIT_ALGEBRAIC_ACSR_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_ACSR_MATRIX_HPP_

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>

namespace NetworKit {

class DenseMatrix;
class VSRMatrix;
template <typename Matrix>
class SMatrixContainer;

/** Represents an augmented matrix entry with @a values.size() values. */
struct ACSRTriplet {
    index row;
    index column;
    std::vector<double> values;
};

/**
 * @ingroup algebraic
 * Sparse matrix stored in augmented compressed sparse row (ACSR) format.
 *
 * ACSR has the same sparsity representation as CSR, but every stored entry has exactly @a k
 * values. The value count is fixed when the matrix is constructed.
 */
class ACSRMatrix final {
    template <typename Matrix>
    friend class SMatrixContainer;

    std::vector<index> rowIdx;
    std::vector<index> columnIdx;
    std::vector<double> nonZeros;

    count nRows;
    count nCols;
    count k;
    bool isSorted;

    static void validateK(count k) {
        if (k == 0)
            throw std::invalid_argument("ACSRMatrix requires k to be positive");
    }

    static void validateValues(const std::vector<double> &values, count k) {
        if (values.size() != k)
            throw std::invalid_argument("Each ACSRMatrix entry must contain exactly k values");
    }

    void assureValues() const {
        if (!columnIdx.empty() && nonZeros.empty())
            throw std::logic_error("This ACSRMatrix stores only a sparsity pattern");
    }

public:
    /** Constructs an empty matrix without an assigned value count. */
    ACSRMatrix()
        : rowIdx(1, 0), columnIdx(), nonZeros(), nRows(0), nCols(0), k(0), isSorted(true) {}

    /** Constructs an empty square matrix whose entries have @a k values. */
    ACSRMatrix(count dimension, count k)
        : rowIdx(dimension + 1, 0), columnIdx(), nonZeros(), nRows(dimension), nCols(dimension),
          k(k), isSorted(true) {
        validateK(k);
    }

    /** Constructs an empty rectangular matrix whose entries have @a k values. */
    ACSRMatrix(count nRows, count nCols, count k)
        : rowIdx(nRows + 1, 0), columnIdx(), nonZeros(), nRows(nRows), nCols(nCols), k(k),
          isSorted(true) {
        validateK(k);
    }

    /** Constructs a square matrix from augmented triplets. */
    ACSRMatrix(count dimension, count k, const std::vector<ACSRTriplet> &triplets)
        : ACSRMatrix(dimension, dimension, k, triplets, false) {}

    /** Constructs a square matrix from augmented triplets. */
    ACSRMatrix(count dimension, const std::vector<ACSRTriplet> &triplets, count k,
               bool isSorted = false)
        : ACSRMatrix(dimension, dimension, k, triplets, isSorted) {}

    /** Constructs a rectangular matrix from augmented triplets. */
    ACSRMatrix(count nRows, count nCols, count k, const std::vector<ACSRTriplet> &triplets,
               bool isSorted = false)
        : rowIdx(nRows + 1, 0), columnIdx(triplets.size()), nonZeros(triplets.size() * k),
          nRows(nRows), nCols(nCols), k(k), isSorted(isSorted) {
        validateK(k);

        for (const auto &triplet : triplets) {
            if (triplet.row >= nRows || triplet.column >= nCols)
                throw std::out_of_range("ACSRMatrix triplet index out of range");
            validateValues(triplet.values, k);
            ++rowIdx[triplet.row + 1];
        }

        for (index i = 0; i < nRows; ++i)
            rowIdx[i + 1] += rowIdx[i];

        auto positions = rowIdx;
        for (const auto &triplet : triplets) {
            const index destination = positions[triplet.row]++;
            columnIdx[destination] = triplet.column;
            std::copy(triplet.values.begin(), triplet.values.end(),
                      nonZeros.begin() + destination * k);
        }
    }

    /** Constructs a rectangular matrix from augmented triplets. */
    ACSRMatrix(count nRows, count nCols, const std::vector<ACSRTriplet> &triplets, count k,
               bool isSorted = false)
        : ACSRMatrix(nRows, nCols, k, triplets, isSorted) {}

    /**
     * Constructs a rectangular matrix from a sparsity pattern and initializes every stored value
     * to @a initialValue. This avoids allocating one value vector per entry when all augmented
     * values are equal.
     */
    ACSRMatrix(count nRows, count nCols, count k, const std::vector<Triplet> &pattern,
               double initialValue, bool isSorted = false)
        : rowIdx(nRows + 1, 0), columnIdx(pattern.size()),
          nonZeros(pattern.size() * k, initialValue), nRows(nRows), nCols(nCols), k(k),
          isSorted(isSorted) {
        validateK(k);

        for (const auto &entry : pattern) {
            if (entry.row >= nRows || entry.column >= nCols)
                throw std::out_of_range("ACSRMatrix pattern index out of range");
            ++rowIdx[entry.row + 1];
        }

        for (index i = 0; i < nRows; ++i)
            rowIdx[i + 1] += rowIdx[i];

        auto positions = rowIdx;
        for (const auto &entry : pattern)
            columnIdx[positions[entry.row]++] = entry.column;
    }

    /** Constructs a matrix from columns and augmented values grouped by row. */
    ACSRMatrix(count nRows, count nCols, count k, const std::vector<std::vector<index>> &columnIdx,
               const std::vector<std::vector<std::vector<double>>> &values, bool isSorted = false)
        : rowIdx(nRows + 1, 0), columnIdx(), nonZeros(), nRows(nRows), nCols(nCols), k(k),
          isSorted(isSorted) {
        validateK(k);
        if (columnIdx.size() != nRows || values.size() != nRows)
            throw std::invalid_argument("ACSRMatrix row data does not match its row count");

        for (index i = 0; i < nRows; ++i) {
            if (columnIdx[i].size() != values[i].size())
                throw std::invalid_argument("ACSRMatrix columns and values have different sizes");

            for (index j = 0; j < columnIdx[i].size(); ++j) {
                if (columnIdx[i][j] >= nCols)
                    throw std::out_of_range("ACSRMatrix column index out of range");
                validateValues(values[i][j], k);
                this->columnIdx.push_back(columnIdx[i][j]);
                nonZeros.insert(nonZeros.end(), values[i][j].begin(), values[i][j].end());
            }
            rowIdx[i + 1] = this->columnIdx.size();
        }
    }

    /** Constructs a matrix from columns and augmented values grouped by row. */
    ACSRMatrix(count nRows, count nCols, const std::vector<std::vector<index>> &columnIdx,
               const std::vector<std::vector<std::vector<double>>> &values, count k,
               bool isSorted = false)
        : ACSRMatrix(nRows, nCols, k, columnIdx, values, isSorted) {}

    /** Constructs a matrix from CSR arrays and one augmented value vector per entry. */
    ACSRMatrix(count nRows, count nCols, count k, const std::vector<index> &rowIdx,
               const std::vector<index> &columnIdx,
               const std::vector<std::vector<double>> &nonZeros, bool isSorted = false)
        : rowIdx(rowIdx), columnIdx(columnIdx), nonZeros(), nRows(nRows), nCols(nCols), k(k),
          isSorted(isSorted) {
        validateK(k);
        if (rowIdx.size() != nRows + 1 || columnIdx.size() != nonZeros.size() || rowIdx.front() != 0
            || rowIdx.back() != nonZeros.size() || !std::is_sorted(rowIdx.begin(), rowIdx.end()))
            throw std::invalid_argument("Invalid ACSRMatrix CSR arrays");

        for (index column : columnIdx) {
            if (column >= nCols)
                throw std::out_of_range("ACSRMatrix column index out of range");
        }
        this->nonZeros.reserve(nonZeros.size() * k);
        for (const auto &values : nonZeros) {
            validateValues(values, k);
            this->nonZeros.insert(this->nonZeros.end(), values.begin(), values.end());
        }
    }

    /** Constructs a matrix from CSR arrays and one augmented value vector per entry. */
    ACSRMatrix(count nRows, count nCols, const std::vector<index> &rowIdx,
               const std::vector<index> &columnIdx,
               const std::vector<std::vector<double>> &nonZeros, count k, bool isSorted = false)
        : ACSRMatrix(nRows, nCols, k, rowIdx, columnIdx, nonZeros, isSorted) {}

    /** Constructs a matrix from CSR arrays and flattened augmented values. */
    ACSRMatrix(count nRows, count nCols, count k, const std::vector<index> &rowIdx,
               const std::vector<index> &columnIdx, const std::vector<double> &nonZeros,
               bool isSorted = false)
        : rowIdx(rowIdx), columnIdx(columnIdx), nonZeros(nonZeros), nRows(nRows), nCols(nCols),
          k(k), isSorted(isSorted) {
        validateK(k);
        if (rowIdx.size() != nRows + 1 || rowIdx.front() != 0 || rowIdx.back() != columnIdx.size()
            || !std::is_sorted(rowIdx.begin(), rowIdx.end())
            || nonZeros.size() != columnIdx.size() * k)
            throw std::invalid_argument("Invalid ACSRMatrix CSR arrays");

        for (index column : columnIdx) {
            if (column >= nCols)
                throw std::out_of_range("ACSRMatrix column index out of range");
        }
    }

    /** Constructs a matrix from CSR arrays and flattened augmented values. */
    ACSRMatrix(count nRows, count nCols, const std::vector<index> &rowIdx,
               const std::vector<index> &columnIdx, const std::vector<double> &nonZeros, count k,
               bool isSorted = false)
        : ACSRMatrix(nRows, nCols, k, rowIdx, columnIdx, nonZeros, isSorted) {}

    /** Constructs a value-bearing matrix by taking ownership of flattened CSR arrays. */
    ACSRMatrix(count nRows, count nCols, count k, std::vector<index> &&rowIdx,
               std::vector<index> &&columnIdx, std::vector<double> &&nonZeros,
               bool isSorted = false)
        : rowIdx(std::move(rowIdx)), columnIdx(std::move(columnIdx)), nonZeros(std::move(nonZeros)),
          nRows(nRows), nCols(nCols), k(k), isSorted(isSorted) {
        validateK(k);
        if (this->rowIdx.size() != nRows + 1 || this->rowIdx.front() != 0
            || this->rowIdx.back() != this->columnIdx.size()
            || !std::is_sorted(this->rowIdx.begin(), this->rowIdx.end())
            || this->nonZeros.size() != this->columnIdx.size() * k)
            throw std::invalid_argument("Invalid ACSRMatrix CSR arrays");

        for (index column : this->columnIdx) {
            if (column >= nCols)
                throw std::out_of_range("ACSRMatrix column index out of range");
        }
    }

    /**
     * Constructs a pattern-only matrix by taking ownership of CSR index arrays. Value-based
     * operations are unavailable, but multiplyInto() can use the pattern directly.
     */
    ACSRMatrix(count nRows, count nCols, count k, std::vector<index> &&rowIdx,
               std::vector<index> &&columnIdx, bool isSorted = false)
        : rowIdx(std::move(rowIdx)), columnIdx(std::move(columnIdx)), nonZeros(), nRows(nRows),
          nCols(nCols), k(k), isSorted(isSorted) {
        validateK(k);
        if (this->rowIdx.size() != nRows + 1 || this->rowIdx.front() != 0
            || this->rowIdx.back() != this->columnIdx.size()
            || !std::is_sorted(this->rowIdx.begin(), this->rowIdx.end()))
            throw std::invalid_argument("Invalid ACSRMatrix pattern arrays");

        for (index column : this->columnIdx) {
            if (column >= nCols)
                throw std::out_of_range("ACSRMatrix column index out of range");
        }
    }

    ACSRMatrix(const ACSRMatrix &other) = default;
    ACSRMatrix(ACSRMatrix &&other) noexcept = default;
    ~ACSRMatrix() = default;

    /**
     * Updates @a other in two passes over the stored entries.
     *
     * The first pass adds each entry's values to the first @a k columns of its source row. The
     * second pass assigns the first @a k columns of each target row to the corresponding entry's
     * values.
     */
    void updateBidirectional(DenseMatrix &other);

    /**
     * Updates @a other in one passes over the stored entries.
     *
     * The first pass adds each entry's values to the first @a k columns of its source row.
     */
    void updateOther(DenseMatrix &other) const;

    void updateOther(VSRMatrix &other) const;

    /**
     * Applies this matrix's sparsity pattern directly to @a input and accumulates the result into
     * @a output. For every stored entry (i,j), the first @a k values of input row j are added to
     * output row i.
     */
    void multiplyInto(const VSRMatrix &input, VSRMatrix &output) const;

    void assign(const DenseMatrix &other);

    void assign(const VSRMatrix &other);

    void resetOther(DenseMatrix &other) const;

    void resetOther(VSRMatrix &other) const;

    ACSRMatrix &operator=(const ACSRMatrix &other) = default;
    ACSRMatrix &operator=(ACSRMatrix &&other) noexcept = default;
};

} // namespace NetworKit

#endif // NETWORKIT_ALGEBRAIC_ACSR_MATRIX_HPP_
