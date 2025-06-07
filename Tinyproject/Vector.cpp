#include "Vector.h"

Vector::Vector(int size) : mSize(size), mData(new double[size]()) {
    assert(size > 0);
}

Vector::Vector(const Vector& other) : mSize(other.mSize), mData(new double[other.mSize]) {
    for (int i = 0; i < mSize; ++i) {
        mData[i] = other.mData[i];
    }
}

Vector::~Vector() {
    delete[] mData;
}

Vector& Vector::operator=(const Vector& other) {
    if (this == &other) return *this;
    if (mSize != other.mSize) {
        delete[] mData;
        mSize = other.mSize;
        mData = new double[mSize];
    }
    for (int i = 0; i < mSize; ++i) {
        mData[i] = other.mData[i];
    }
    return *this;
}

Vector Vector::operator-() const {
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result.mData[i] = -mData[i];
    }
    return result;
}

Vector Vector::operator+(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result.mData[i] = mData[i] + other.mData[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result.mData[i] = mData[i] - other.mData[i];
    }
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result.mData[i] = mData[i] * scalar;
    }
    return result;
}

double& Vector::operator[](int index) {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

const double& Vector::operator[](int index) const {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

double& Vector::operator()(int index) {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

const double& Vector::operator()(int index) const {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

int Vector::GetSize() const {
    return mSize;
}

Vector operator*(double scalar, const Vector& v) {
    return v * scalar;
}

std::ostream& operator<<(std::ostream& os, const Vector& v) {
    os << std::fixed << std::setprecision(6) << "[";
    for (int i = 0; i < v.mSize; ++i) {
        os << v.mData[i];
        if (i < v.mSize - 1) os << ", ";
    }
    os << "]";
    return os;
}

std::istream& operator>>(std::istream& is, Vector& v) {
    for (int i = 0; i < v.mSize; ++i) {
        is >> v.mData[i];
    }
    return is;
}

