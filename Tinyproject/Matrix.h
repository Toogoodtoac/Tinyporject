#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

class Matrix {
private:
    int    mNumRows;
    int    mNumCols;
    double** mData;

public:
    Matrix(int numRows, int numCols);
    Matrix(const Matrix& other);
    ~Matrix();

    Matrix& operator=(const Matrix& other);
    Matrix  operator-() const;
    Matrix  operator+(const Matrix& other) const;
    Matrix  operator-(const Matrix& other) const;
    Matrix  operator*(const Matrix& other) const;
    Matrix  operator*(double scalar) const;
    Vector  operator*(const Vector& v) const;

    double&       operator()(int row, int col);
    const double& operator()(int row, int col) const;

    int GetNumRows() const;
    int GetNumCols() const;

    double computeDeterminant() const;
    Matrix computeInverse() const;
    Matrix computePseudoInverse() const;
    Matrix Transpose() const;

    friend Matrix        operator*(double scalar, const Matrix& m);
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    friend std::istream& operator>>(std::istream& is, Matrix& m);
};

#endif // MATRIX_H
