#ifndef POSSYMLINSYSTEM_H
#define POSSYMLINSYSTEM_H

#include "LinearSystem.h"
#include <cassert>
#include <cmath>

class PosSymLinSystem : public LinearSystem {
public:
    PosSymLinSystem(const Matrix& A, const Vector& b);
    Vector Solve() override;

private:
    bool IsSymmetric(const Matrix& A) const;
};

#endif // POSSYMLINSYSTEM_H

