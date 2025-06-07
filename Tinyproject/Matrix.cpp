#include "Matrix.h"

// Constructor: allocate & zero‐initialize
Matrix::Matrix(int numRows, int numCols)
    : mNumRows(numRows), mNumCols(numCols)
{
    assert(numRows > 0 && numCols > 0);
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; ++i) {
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; ++j) {
            mData[i][j] = 0.0;
        }
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& otherMatrix)
    : mNumRows(otherMatrix.mNumRows), mNumCols(otherMatrix.mNumCols)
{
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; ++i) {
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; ++j) {
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }
}

// Destructor
Matrix::~Matrix() {
    for (int i = 0; i < mNumRows; ++i) {
        delete[] mData[i];
    }
    delete[] mData;
}

// Assignment operator
Matrix& Matrix::operator=(const Matrix& otherMatrix) {
    if (this == &otherMatrix) return *this;

    // If dimensions differ, reallocate
    if (mNumRows != otherMatrix.mNumRows || mNumCols != otherMatrix.mNumCols) {
        for (int i = 0; i < mNumRows; ++i) {
            delete[] mData[i];
        }
        delete[] mData;

        mNumRows = otherMatrix.mNumRows;
        mNumCols = otherMatrix.mNumCols;
        mData = new double*[mNumRows];
        for (int i = 0; i < mNumRows; ++i) {
            mData[i] = new double[mNumCols];
        }
    }
    // Copy elements
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            mData[i][j] = otherMatrix.mData[i][j];
        }
    }
    return *this;
}

// Unary negation (element‐wise)
Matrix Matrix::operator-() const {
    Matrix m(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            m.mData[i][j] = -mData[i][j];
        }
    }
    return m;
}

// Matrix + Matrix
Matrix Matrix::operator+(const Matrix& otherMatrix) const {
    assert(mNumRows == otherMatrix.mNumRows && mNumCols == otherMatrix.mNumCols);
    Matrix m(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            m.mData[i][j] = mData[i][j] + otherMatrix.mData[i][j];
        }
    }
    return m;
}

// Matrix - Matrix
Matrix Matrix::operator-(const Matrix& otherMatrix) const {
    assert(mNumRows == otherMatrix.mNumRows && mNumCols == otherMatrix.mNumCols);
    Matrix m(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            m.mData[i][j] = mData[i][j] - otherMatrix.mData[i][j];
        }
    }
    return m;
}

// Matrix * Matrix
Matrix Matrix::operator*(const Matrix& otherMatrix) const {
    assert(mNumCols == otherMatrix.mNumRows);
    Matrix m(mNumRows, otherMatrix.mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < otherMatrix.mNumCols; ++j) {
            double sum = 0.0;
            for (int k = 0; k < mNumCols; ++k) {
                sum += mData[i][k] * otherMatrix.mData[k][j];
            }
            m.mData[i][j] = sum;
        }
    }
    return m;
}

// Matrix * Vector
Vector Matrix::operator*(const Vector& v) const {
    assert(mNumCols == v.GetSize());
    Vector result(mNumRows);
    for (int i = 0; i < mNumRows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < mNumCols; ++j) {
            sum += mData[i][j] * v[j]; // zero‐based on v
        }
        result[i] = sum; // zero‐based on result
    }
    return result;
}

// Matrix * scalar
Matrix Matrix::operator*(double scalar) const {
    Matrix m(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            m.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return m;
}

// One‐based indexing (read/write)
double& Matrix::operator()(int row, int col) {
    assert(row >= 1 && row <= mNumRows && col >= 1 && col <= mNumCols);
    return mData[row - 1][col - 1];
}

// One‐based indexing (read‐only)
const double& Matrix::operator()(int row, int col) const {
    assert(row >= 1 && row <= mNumRows && col >= 1 && col <= mNumCols);
    return mData[row - 1][col - 1];
}

// Accessors
int Matrix::GetNumRows() const { return mNumRows; }
int Matrix::GetNumCols() const { return mNumCols; }

// Determinant (Gaussian elimination approach)
double Matrix::computeDeterminant() const {
    assert(mNumRows == mNumCols);
    if (mNumRows == 1) {
        return mData[0][0];
    }

    Matrix tempMatrix = *this; // make a copy
    double det = 1.0;
    int numSwaps = 0;

    for (int k = 0; k < mNumRows; ++k) {
        // find pivot row (max absolute value in column k at or below row k)
        int pivotRow = k;
        for (int r = k + 1; r < mNumRows; ++r) {
            if (std::fabs(tempMatrix.mData[r][k]) > std::fabs(tempMatrix.mData[pivotRow][k])) {
                pivotRow = r;
            }
        }
        // if pivotRow ≠ k, swap entire row
        if (pivotRow != k) {
            std::swap(tempMatrix.mData[k], tempMatrix.mData[pivotRow]);
            numSwaps++;
        }
        // if pivot is (nearly) zero → determinant = 0
        if (std::fabs(tempMatrix.mData[k][k]) < DOUBLE_TOLERANCE) {
            return 0.0;
        }
        det *= tempMatrix.mData[k][k];

        // eliminate below
        for (int i = k + 1; i < mNumRows; ++i) {
            double factor = tempMatrix.mData[i][k] / tempMatrix.mData[k][k];
            for (int j = k; j < mNumCols; ++j) {
                tempMatrix.mData[i][j] -= factor * tempMatrix.mData[k][j];
            }
        }
    }
    if (numSwaps % 2 != 0) {
        det = -det;
    }
    return det;
}

// Inverse via Gaussian elimination (augmented matrix [A | I])
Matrix Matrix::computeInverse() const {
    assert(mNumRows == mNumCols);
    int n = mNumRows;
    Matrix augmented(n, 2 * n);

    // build [A | I]
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented.mData[i][j] = mData[i][j];
        }
        augmented.mData[i][i + n] = 1.0;
    }

    // forward elimination with pivoting
    for (int k = 0; k < n; ++k) {
        // find pivot in column k (largest absolute value)
        int pivotRow = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::fabs(augmented.mData[i][k]) > std::fabs(augmented.mData[pivotRow][k])) {
                pivotRow = i;
            }
        }
        // swap if needed
        if (pivotRow != k) {
            std::swap(augmented.mData[k], augmented.mData[pivotRow]);
        }
        // if pivot is zero → singular
        if (std::fabs(augmented.mData[k][k]) < DOUBLE_TOLERANCE) {
            std::cerr << "Error: Matrix is singular or ill-conditioned; cannot invert.\n";
            return Matrix(n, n); // return a zero matrix (you could also throw)
        }

        // normalize pivot row
        double pivotVal = augmented.mData[k][k];
        for (int j = k; j < 2 * n; ++j) {
            augmented.mData[k][j] /= pivotVal;
        }
        // eliminate all other rows
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            double factor = augmented.mData[i][k];
            for (int j = k; j < 2 * n; ++j) {
                augmented.mData[i][j] -= factor * augmented.mData[k][j];
            }
        }
    }

    // extract inverse from the right half
    Matrix inverse(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse.mData[i][j] = augmented.mData[i][j + n];
        }
    }
    return inverse;
}

// Pseudo‐inverse via normal equations: A⁺ = (Aᵀ A)⁻¹ Aᵀ
Matrix Matrix::computePseudoInverse() const {
    Matrix A_T = this->Transpose();   // (numCols × numRows)
    Matrix A_T_A = A_T * (*this);   // (numCols × numRows) * (numRows × numCols) = (numCols × numCols)
    Matrix A_T_A_inv = A_T_A.computeInverse(); // invert
    return A_T_A_inv * A_T;    // (numCols × numCols) * (numCols × numRows) = (numCols × numRows)
}

// Transpose
Matrix Matrix::Transpose() const {
    Matrix transposed(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            transposed.mData[j][i] = mData[i][j];
        }
    }
    return transposed;
}

// Friend: scalar * Matrix
Matrix operator*(double scalar, const Matrix& m) {
    return m * scalar;
}

//  print Matrix
std::ostream& operator<<(std::ostream& output, const Matrix& m) {
    output << std::fixed << std::setprecision(6);
    for (int i = 0; i < m.mNumRows; ++i) {
        output << "[";
        for (int j = 0; j < m.mNumCols; ++j) {
            output << std::setw(10) << m.mData[i][j];
            if (j < m.mNumCols - 1) {
                output << ", ";
            }
        }
        output << "]\n";
    }
    return output;
}

// read Matrix elements
std::istream& operator>>(std::istream& input, Matrix& m) {
    std::cout << "Enter " << m.mNumRows << "x" << m.mNumCols << " elements:\n";
    for (int i = 0; i < m.mNumRows; ++i) {
        for (int j = 0; j < m.mNumCols; ++j) {
            input >> m.mData[i][j];
        }
    }
    return input;
}

