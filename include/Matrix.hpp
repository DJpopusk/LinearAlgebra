#pragma once
#include <vector>
#include <thread>
#include "Exceptions.hpp"

class Matrix {
public:
    Matrix(size_t rows, size_t cols);
    Matrix(const Matrix& other) = default;
    Matrix& operator=(const Matrix& other) = default;

    size_t getRows() const { return rows_; }
    size_t getCols() const { return cols_; }

    double& operator()(size_t row, size_t col);
    const double& operator()(size_t row, size_t col) const;

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;

    friend Matrix operator*(double scalar, const Matrix& matrix);

private:
    void distribute_rows(size_t total_rows, size_t num_threads,
                        std::vector<std::thread>& threads,
                        std::function<void(size_t, size_t)> worker) const;
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};