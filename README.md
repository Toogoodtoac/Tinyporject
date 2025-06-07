# Linear Regression Project

This project implements core linear algebra structures in C++ and applies them to a real-world regression problem: predicting CPU performance from hardware features.

## 🔧 What Was Done
- Built custom `Vector`, `Matrix`, and `LinearSystem` classes from scratch.
- Implemented solvers:
  - Gaussian Elimination
  - Conjugate Gradient (for symmetric positive-definite systems)
  - Moore-Penrose Pseudo-inverse
  - Tikhonov Regularization
- Used these solvers to perform linear regression on the [UCI CPU Performance dataset](https://archive.ics.uci.edu/ml/datasets/Computer+Hardware).

## 📊 Problem
Predict Published Relative Performance (PRP) from 6 hardware features:
- MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX

## 🧠 Approach
- Removed non-numeric columns.
- Split dataset into 80% train, 20% test.
- Solved using normal equation and Tikhonov regularized variant.

## 📈 Results
- **Moore-Penrose Pseudo-inverse RMSE**: ~59.86
- **Tikhonov Regularization RMSE**: ~59.86
- Performance was limited due to model simplicity and dataset characteristics.

## 📁 Files
- `Vector.hpp`, `Matrix.hpp`, `LinearSystem.hpp`: Custom data structures
- `main.cpp`: Training + prediction pipeline
- `report.tex`: Full LaTeX project report

## ✅ To Run
```bash
g++ -std=c++17 main.cpp -o linear_regression
./linear_regression
