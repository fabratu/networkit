#include <algorithm>
#include <cassert>

#include <networkit/algebraic/ACSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

void ACSRMatrix::updateOther(DenseMatrix &other) const {
    for (index i = 0; i < nRows; ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index valueOffset = entry * k;
            for (index value = 0; value < k; ++value)
                other(i, value) += nonZeros[valueOffset + value];
        }
    }
}

void ACSRMatrix::assign(const DenseMatrix &other) {
    for (index i = 0; i < nRows; ++i) {
        for (index entry = rowIdx[i]; entry < rowIdx[i + 1]; ++entry) {
            const index j = columnIdx[entry];
            const index valueOffset = entry * k;
            for (index value = 0; value < k; ++value)
                nonZeros[valueOffset + value] = other(j, value);
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

void ACSRMatrix::updateBidirectional(DenseMatrix &other) {
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
