# Tiny Project 

## Overview

This project implements custom C++ classes for solving systems of linear equations and applies them to predict CPU performance using linear regression. It includes support for both square and non-square systems, as well as methods suitable for symmetric positive-definite matrices.

## Structure

### Part A: Linear Solvers

- `Vector.cpp`: Implements a 1D vector with operator overloading and indexing.
- `Matrix.cpp`: Implements a 2D matrix with support for basic operations, determinant, inverse, and pseudo-inverse.
- `LinearSystem.cpp`: Solves square systems using Gaussian elimination with pivoting.
- `PosSymLinSystem.cpp`: Solves symmetric positive-definite systems using the **Conjugate Gradient method**.

### Part B: Regression and Real-world Application

- `main.cpp`: Applies **Moore-Penrose Pseudo-inverse** and **Tikhonov Regularization** for linear regression.
    - Loads and processes CPU performance dataset.
    - Trains regression models using your solvers.
    - Evaluates performance using RMSE.

## Results

- RMSE (Pseudo-inverse): ~59.86
- RMSE (Tikhonov): ~59.86
- Conjugate Gradient showed fast convergence on SPD systems

## Author

Ngo Anh Hieu — Vietnamese-German University
