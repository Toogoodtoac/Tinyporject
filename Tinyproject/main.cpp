#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <iomanip> // For std::fixed and std::setprecision
#include <chrono>

#include "Vector.h"
#include "Matrix.h"
#include "LinearSystem.h"
// #include "PosSymLinSystem.h" // Uncomment to solve SPD systems with CG

// Function to read data from the machine.data CSV file
// Skips the first two non-predictive attributes (vendor name, model name)
// Extracts MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX for features (X) and PRP for target (y)
bool read_data(const std::string& filename,
               std::vector<std::vector<double>>& X_data,
               std::vector<double>& y_data)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        int         idx = 0;
        std::vector<double> features;
        double target = 0.0;

        while (std::getline(ss, token, ',')) {
            token.erase(0, token.find_first_not_of(" \t\n\r\f\v"));
            token.erase(token.find_last_not_of(" \t\n\r\f\v") + 1);

            if (idx >= 2) {
                try {
                    double val = std::stod(token);
                    if
                    (
                        idx == 2 || // MYCT
                        idx == 3 || // MMIN
                        idx == 4 || // MMAX
                        idx == 5 || // CACH
                        idx == 6 || // CHMIN <--- ADDED CHMIN
                        idx == 7)   // CHMAX
                    {
                        features.push_back(val);
                    } else if (idx == 8)
                    { // PRP (target variable)
                        target = val;
                    }
                }
                catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid argument: " << token << " in line: " << line << " - " << e.what() << "\n";
                    return false;
                } catch (const std::out_of_range& e) {
                    std::cerr << "Out of range: " << token << " in line: " << line << " - " << e.what() << "\n";
                    return false;
                }
            }
            ++idx;
        }
        if (features.size() == 6) { // Ensure all 6 features were collected <--- UPDATED FROM 5 TO 6
            X_data.push_back(features);
            y_data.push_back(target);
        }
    }
    file.close();
    return true;
}

// Function to split data into training and testing sets (80/20 split)
void split_data(
    const std::vector<std::vector<double>>& X_data,
    const std::vector<double>& y_data,
    std::vector<std::vector<double>>& X_train,
    std::vector<double>& y_train,
    std::vector<std::vector<double>>& X_test,
    std::vector<double>& y_test,
    double train_ratio = 0.8)
{
    int N = (int)X_data.size();
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, 2, ... N-1

    // Randomly shuffle the indices for reproducible split
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed); // Mersenne Twister engine
    std::shuffle(indices.begin(), indices.end(), rng);

    int train_size = static_cast<int>(N * train_ratio);
    for (int i = 0; i < train_size; ++i) {
        X_train.push_back(X_data[indices[i]]);
        y_train.push_back(y_data[indices[i]]);
    }
    for (int i = train_size; i < N; ++i) {
        X_test.push_back(X_data[indices[i]]);
        y_test.push_back(y_data[indices[i]]);
    }
}

// Function to calculate Root Mean Square Error (RMSE)
double calculate_rmse(const std::vector<double>& actual,const std::vector<double>& predicted)
{
    assert(actual.size() == predicted.size());
    if (actual.empty()) return 0.0;
    double sum_sq = 0.0;
    for (size_t i = 0; i < actual.size(); ++i) {
        double diff = actual[i] - predicted[i];
        sum_sq += diff * diff;
    }
    return std::sqrt(sum_sq / actual.size());
}

// Function to solve linear regression using Tikhonov Regularization
// Solves (A^T A + mu * I) x = A^T b
Vector solve_tikhonov_regularization(const Matrix& A, const Vector& b, double mu) {
    assert(mu >= 0.0); // Regularization parameter must be non-negative

    Matrix A_T = A.Transpose();
    Matrix A_T_A = A_T * A;
    Vector A_T_b = A_T * b;

    int n_features = A.GetNumCols();
    // Create an identity matrix I of size n_features x n_features
    Matrix I(n_features, n_features);
    for (int i = 1; i <= n_features; ++i) {
        I(i, i) = 1.0;
    }

    // Form the regularized matrix (A^T A + mu * I)
    Matrix regularized_matrix = A_T_A + (I * mu);

    // Solve the regularized system (regularized_matrix * x = A_T_b)
    LinearSystem system(regularized_matrix, A_T_b);
    Vector x = system.Solve();
    return x;
}


int main() {
    std::cout << "--- CPU Performance Prediction (Linear Regression) ---\n";

    std::vector<std::vector<double>> X_raw; // Features: MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX
    std::vector<double> y_raw;              // Target: PRP
    std::string filename = "machine.data";

    // 1. Data Acquisition and Preprocessing
    if (!read_data(filename, X_raw, y_raw)) {
        std::cerr << "Failed to read data from " << filename << ". Please ensure the file exists." << "\n";
        return 1;
    }
    std::cout << "Loaded " << X_raw.size() << " instances from " << filename << "\n";

    std::vector<std::vector<double>> X_train, X_test;
    std::vector<double> y_train, y_test;
    split_data(X_raw, y_raw, X_train, y_train, X_test, y_test, 0.8);
    std::cout << "Data split: " << X_train.size() << " training instances, "<< X_test.size() << " testing instances." << "\n";

    // 2. Constructing the Design Matrix (A) and Target Vector (b)
    int n_train = (int)X_train.size();
    const int features = 6; // MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX
    Matrix A_train(n_train, features);
    Vector b_train(n_train);

    for (int i = 0; i < n_train; ++i) {
        // Populate A_train with features (using 1-based indexing for Matrix operator())
        A_train(i + 1, 1) = X_train[i][0]; // MYCT
        A_train(i + 1, 2) = X_train[i][1]; // MMIN
        A_train(i + 1, 3) = X_train[i][2]; // MMAX
        A_train(i + 1, 4) = X_train[i][3]; // CACH
        A_train(i + 1, 5) = X_train[i][4]; // CHMIN
        A_train(i + 1, 6) = X_train[i][5]; // CHMAX

        // Populate b_train with target PRP values (using 1-based indexing for Vector operator())
        b_train(i + 1) = y_train[i];
    }
    std::cout << "\nConstructed design matrix A (" << A_train.GetNumRows() << "x" << A_train.GetNumCols()<< ") and target vector b (" << b_train.GetSize() << ").\n";


    // 3. Solving the Linear System for Parameters (x) using Moore-Penrose Pseudo-inverse
    std::cout << "\n--- Solving using Moore-Penrose Pseudo-inverse ---" << "\n";
    std::cout << "Computing pseudo-inverse of A_train …\n";
    Matrix A_pinv = A_train.computePseudoInverse();
    std::cout << "Pseudo-inverse computation complete." << "\n";

    Vector x_params_pinv = A_pinv * b_train;
    std::cout << "Learned Parameters (Pseudo-inverse): " << x_params_pinv << "\n";

    // 4. Model Evaluation on the Testing Set (Pseudo-inverse)
    std::vector<double> predicted_test_pinv;
    for (size_t i = 0; i < X_test.size(); ++i) {
        double y_hat = 0.0;
        for (int j = 0; j < features; ++j) { // Loop for all 6 features
            y_hat += x_params_pinv[j] * X_test[i][j];
        }
        predicted_test_pinv.push_back(y_hat);
    }
    double rmse_pinv = calculate_rmse(y_test, predicted_test_pinv);
    std::cout << "RMSE on test set (Pseudo-inverse): "
              << std::fixed << std::setprecision(6) << rmse_pinv << "\n";


    // 5. Solving the Linear System for Parameters (x) using Tikhonov Regularization
    std::cout << "\n--- Solving using Tikhonov Regularization ---" << "\n";
    double mu = 0.1; // Example regularization parameter
    std::cout << "Using regularization parameter (mu) = " << mu << "\n";

    Vector x_params_tikhonov = solve_tikhonov_regularization(A_train, b_train, mu);
    std::cout << "Learned Parameters (Tikhonov Regularization): " << x_params_tikhonov << "\n";

    // 6. Model Evaluation on the Testing Set (Tikhonov Regularization)
    std::vector<double> predicted_test_tikhonov;
    for (size_t i = 0; i < X_test.size(); ++i) {
        double y_hat = 0.0;
        for (int j = 0; j < features; ++j) { // Loop for all 6 features
            y_hat += x_params_tikhonov[j] * X_test[i][j];
        }
        predicted_test_tikhonov.push_back(y_hat);
    }
    double rmse_tikhonov = calculate_rmse(y_test, predicted_test_tikhonov);
    std::cout << "RMSE on test set (Tikhonov Regularization): "<< std::fixed << std::setprecision(6) << rmse_tikhonov << "\n";

    std::cout << "\n--- Program Finished ---" << "\n";
    return 0;
}
