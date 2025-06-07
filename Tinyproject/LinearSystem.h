#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Matrix.h"
#include "Vector.h"
#include <cassert>
#include <cmath>

class LinearSystem {
protected:
    int    mSize;
    Matrix* mpA;
    Vector* mpb;

public:
    explicit    LinearSystem(const Matrix& A, const Vector& b);
    virtual    ~LinearSystem();

    virtual Vector Solve();
    int        GetSize() const;
};

#endif // LINEARSYSTEM_H

