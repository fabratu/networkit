#ifndef NETWORKIT_ALGEBRAIC_VSR_MATRIX_HPP_
#define NETWORKIT_ALGEBRAIC_VSR_MATRIX_HPP_

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/Vector.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Rectangular variable sparse row matrix.
 *
 * Row @c i stores a dense prefix of @c kValues[i] entries starting at column zero. Every prefix
 * length is between one and @a k, and @a k must be smaller than the number of columns. Stored
 * values are initialized to 0.0.
 */
class VSRMatrix final {
    std::vector<index> rowIdx;
    std::vector<double> values;

    count nRows{0};
    count nCols{0};
    count k{0};

    static count inferK(const std::vector<count> &kValues);

public:
    /** Constructs an empty matrix. */
    VSRMatrix() = default;

    /** Constructs a matrix with the supplied prefix length for every row. */
    VSRMatrix(count nRows, count nCols, count k, const std::vector<count> &kValues);

    /** Constructs a matrix with @a k after the row-prefix lengths. */
    VSRMatrix(count nRows, count nCols, const std::vector<count> &kValues, count k)
        : VSRMatrix(nRows, nCols, k, kValues) {}

    /** Constructs a matrix and infers @a k as the largest supplied prefix length. */
    VSRMatrix(count nRows, count nCols, const std::vector<count> &kValues)
        : VSRMatrix(nRows, nCols, inferK(kValues), kValues) {}

    /** Returns the value at position (@a i, @a j), or zero outside the stored row prefix. */
    const double &operator()(index i, index j) const;

    /** Returns a mutable reference to a value in the stored prefix of row @a i. */
    double &operator()(index i, index j);

    /** Resets all entries to zero */
    void reset();

    /** Multiplies this matrix with a column vector. */
    Vector operator*(const Vector &vector) const;

    /** Multiplies this matrix with a column vector in-place */
    void muVInPlace(const Vector &vector, Vector &other) const;
};

} // namespace NetworKit

#endif // NETWORKIT_ALGEBRAIC_VSR_MATRIX_HPP_
