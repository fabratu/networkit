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

    values.resize(rowIdx.back(), 0.0);
}

Vector VSRMatrix::operator*(const Vector &vector) const {
    assert(!vector.isTransposed());
    assert(vector.getDimension() == nCols);

    Vector result(nRows, 0.0);
#pragma omp parallel for
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
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(nRows); ++i) {
        double sum = 0.0;
        const index rowBegin = rowIdx[i];
        const index rowEnd = rowIdx[i + 1];
        for (index entry = rowBegin; entry < rowEnd; ++entry)
            sum += values[entry] * vector[entry - rowBegin];
        other[i] = sum;
    }
}

} // namespace NetworKit
