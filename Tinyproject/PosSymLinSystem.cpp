
#include "PosSymLinSystem.h"
#include <iostream>

PosSymLinSystem::PosSymLinSystem(const Matrix& A, const Vector& b)
    : LinearSystem(A, b)
{
    assert(IsSymmetric(A));
}

bool PosSymLinSystem::IsSymmetric(const Matrix& A) const {
    int n = A.GetNumRows();
    if (n != A.GetNumCols()) return false;
    for (int i = 1; i <= n; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            if (std::fabs(A(i, j) - A(j, i)) > DOUBLE_TOLERANCE) {
                return false;
            }
        }
    }
    return true;
}

Vector PosSymLinSystem::Solve() {
    const Matrix& A = *mpA;
    const Vector& b = *mpb;
    int n = mSize;

    Vector x(n);
    Vector r = b;
    Vector p = r;

    double r_dot_r = 0.0;
    for (int i = 0; i < n; ++i) {
        r_dot_r += r[i] * r[i];
    }
    if (std::sqrt(r_dot_r) < DOUBLE_TOLERANCE) {
        return x;
    }

    double tol = 1e-6;
    int maxIt = 2 * n;

    for (int k = 0; k < maxIt; ++k) {
        Vector Ap = A * p;
        double p_dot_Ap = 0.0;
        for (int i = 0; i < n; ++i) {
            p_dot_Ap += p[i] * Ap[i];
        }
        if (std::fabs(p_dot_Ap) < DOUBLE_TOLERANCE) {
            std::cerr << "Warning: CG denominator nearly zero.\n";
            break;
        }

        double alpha = r_dot_r / p_dot_Ap;
        x = x + (p * alpha);
        r = r - (Ap * alpha);

        double r_dot_r_new = 0.0;
        for (int i = 0; i < n; ++i) {
            r_dot_r_new += r[i] * r[i];
        }
        if (std::sqrt(r_dot_r_new) < tol) {
            break;
        }

        double beta = r_dot_r_new / r_dot_r;
        p = r + (p * beta);
        r_dot_r = r_dot_r_new;
    }
    return x;
}
