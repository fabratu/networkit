#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <networkit/algebraic/VSRMatrix.hpp>

namespace NetworKit {

count VSRMatrix::inferK(const std::vector<count> &kValues) {
    if (kValues.empty())
        throw std::invalid_argument("VSRMatrix cannot infer k from an empty row-length vector");
    return *std::max_element(kValues.begin(), kValues.end());
}

VSRMatrix::VSRMatrix(count nRows, count nCols, count k, const std::vector<count> &kValues)
    : rowIdx(nRows + 1, 0), values(), nRows(nRows), nCols(nCols), k(k) {
    if (kValues.size() != nRows)
        throw std::invalid_argument("VSRMatrix requires one prefix length per row");
    if (k == 0 || k >= nCols)
        throw std::invalid_argument("VSRMatrix requires 0 < k < number of columns");

    for (index i = 0; i < nRows; ++i) {
        if (kValues[i] == 0 || kValues[i] > k)
            throw std::invalid_argument("VSRMatrix row prefix lengths must be between 1 and k");
        rowIdx[i + 1] = rowIdx[i] + kValues[i];
    }

    values.resize(rowIdx.back(), 1.0);
}

const double &VSRMatrix::operator()(index i, index j) const {
    if (i >= nRows || j >= nCols)
        throw std::out_of_range("VSRMatrix index out of range");

    static const double zero = 0.0;
    const index rowLength = rowIdx[i + 1] - rowIdx[i];
    return j < rowLength ? values[rowIdx[i] + j] : zero;
}

double &VSRMatrix::operator()(index i, index j) {
    if (i >= nRows || j >= nCols)
        throw std::out_of_range("VSRMatrix index out of range");
    if (j >= rowIdx[i + 1] - rowIdx[i])
        throw std::out_of_range("VSRMatrix index is outside the stored row prefix");

    return values[rowIdx[i] + j];
}

Vector VSRMatrix::operator*(const Vector &vector) const {
    assert(!vector.isTransposed());
    assert(vector.getDimension() == nCols);

    Vector result(nRows, 0.0);
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        double sum = 0.0;
        const index rowBegin = rowIdx[i];
        const index rowEnd = rowIdx[i + 1];
        for (index entry = rowBegin; entry < rowEnd; ++entry)
            sum += values[entry] * vector[entry - rowBegin];
        result[i] = sum;
    }

    return result;
}

void VSRMatrix::muVInPlace(const Vector &vector, Vector &other) const {
#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        double sum = 0.0;
        const index rowBegin = rowIdx[i];
        const index rowEnd = rowIdx[i + 1];
        for (index entry = rowBegin; entry < rowEnd; ++entry)
            sum += values[entry] * vector[entry - rowBegin];
        other[i] = sum;
    }
}

void VSRMatrix::reset() {
    std::fill(values.begin(), values.end(), 0.0);
}

} // namespace NetworKit
