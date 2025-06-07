
#include "LinearSystem.h"
#include <iostream>

LinearSystem::LinearSystem(const Matrix& A, const Vector& b)
    : mSize(A.GetNumRows()), mpA(new Matrix(A)), mpb(new Vector(b))
{
    assert(A.GetNumRows() == A.GetNumCols());
    assert(A.GetNumRows() == b.GetSize());
}

LinearSystem::~LinearSystem() {
    delete mpA;
    delete mpb;
}

Vector LinearSystem::Solve() {
    Matrix A_local = *mpA;
    Vector b_local = *mpb;
    int n = mSize;

    for (int k = 0; k < n; ++k) {
        int pivot = k;
        for (int r = k + 1; r < n; ++r) {
            if (std::fabs(A_local(r + 1, k + 1)) > std::fabs(A_local(pivot + 1, k + 1))) {
                pivot = r;
            }
        }
        if (pivot != k) {
            for (int c = 0; c < n; ++c) {
                std::swap(A_local(k + 1, c + 1), A_local(pivot + 1, c + 1));
            }
            std::swap(b_local(k + 1), b_local(pivot + 1));
        }
        if (std::fabs(A_local(k + 1, k + 1)) < DOUBLE_TOLERANCE) {
            std::cerr << "Error: singular or ill-conditioned matrix.\n";
            return Vector(n);
        }
        for (int r = k + 1; r < n; ++r) {
            double factor = A_local(r + 1, k + 1) / A_local(k + 1, k + 1);
            for (int c = k; c < n; ++c) {
                A_local(r + 1, c + 1) -= factor * A_local(k + 1, c + 1);
            }
            b_local(r + 1) -= factor * b_local(k + 1);
        }
    }

    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A_local(i + 1, j + 1) * x(j + 1);
        }
        x(i + 1) = (b_local(i + 1) - sum) / A_local(i + 1, i + 1);
    }
    return x;
}

int LinearSystem::GetSize() const {
    return mSize;
}
