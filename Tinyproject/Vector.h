
#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cassert>
#include <iomanip>

const double DOUBLE_TOLERANCE = 1e-9;

class Vector {
private:
    int    mSize;
    double* mData;

public:
    explicit Vector(int size);
    Vector(const Vector& other);
    ~Vector();

    Vector& operator=(const Vector& other);
    Vector  operator-() const;
    Vector  operator+(const Vector& other) const;
    Vector  operator-(const Vector& other) const;
    Vector  operator*(double scalar) const;

    double& operator[](int index);
    const double& operator[](int index) const;
    double& operator()(int index);
    const double& operator()(int index) const;

    int GetSize() const;

    friend Vector operator*(double scalar, const Vector& v);
    friend std::ostream& operator<<(std::ostream& os, const Vector& v);
    friend std::istream& operator>>(std::istream& is, Vector& v);
};

#endif // VECTOR_H
