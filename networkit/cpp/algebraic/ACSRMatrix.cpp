#include <algorithm>
#include <cassert>
#include <functional>

#include <networkit/algebraic/ACSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/VSRMatrix.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

void ACSRMatrix::updateOther(DenseMatrix &other) const {
    assureValues();
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index valueOffset = entry * k;
            for (index value = 0; value < k; ++value)
                other(i, value) += nonZeros[valueOffset + value];
        }
    }
}

void ACSRMatrix::updateOther(VSRMatrix &other) const {
    assureValues();
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index valueOffset = entry * k;
            double *const otherRow = &other(i, 0);
            std::transform(nonZeros.cbegin() + valueOffset, nonZeros.cbegin() + valueOffset + k,
                           otherRow, otherRow, std::plus<double>{});
        }
    }
}

void ACSRMatrix::multiplyInto(const VSRMatrix &input, VSRMatrix &output) const {
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        if (rowIdx[i] == rowIdx[i + 1])
            continue;
        double *const outputRow = &output(i, 0);
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const double *const inputRow = &input(columnIdx[entry], 0);
            std::transform(inputRow, inputRow + k, outputRow, outputRow, std::plus<double>{});
        }
    }
}

void ACSRMatrix::assign(const DenseMatrix &other) {
    assureValues();
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index j = columnIdx[entry];
            const index valueOffset = entry * k;
            for (index value = 0; value < k; ++value)
                nonZeros[valueOffset + value] = other(j, value);
        }
    }
}

void ACSRMatrix::assign(const VSRMatrix &other) {
    assureValues();
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index j = columnIdx[entry];
            const index valueOffset = entry * k;
            std::copy_n(&other(j, 0), k, nonZeros.begin() + valueOffset);
        }
    }
}

void ACSRMatrix::resetOther(DenseMatrix &other) const {
    for (index i = 0; i < nRows; ++i) {
        for (index value = 0; value < k; ++value) {
            other(i, value) = 0.0;
        }
    }
}

void ACSRMatrix::resetOther(VSRMatrix &other) const {
    other.reset();
}

void ACSRMatrix::updateBidirectional(DenseMatrix &other) {
    assureValues();
    assert(other.numberOfRows() >= std::max(nRows, nCols));
    assert(other.numberOfColumns() >= k);

    for (index i = 0; i < nRows; ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index valueOffset = entry * k;
            for (index value = 0; value < k; ++value)
                other(i, value) += nonZeros[valueOffset + value];
        }
    }

    for (index i = 0; i < nRows; ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index j = columnIdx[entry];
            const index valueOffset = entry * k;
            for (index value = 0; value < k; ++value)
                nonZeros[valueOffset + value] = other(j, value);
        }
    }
}

} // namespace NetworKit
