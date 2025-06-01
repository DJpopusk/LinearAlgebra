#pragma once
#include <cstddef>
#include <iostream>
#include <thread>
#include <vector>
#include <functional>
#include "Exceptions.hpp"
#include "ThreadPoolConfig.hpp"

class Vector;

class Matrix {
public:
    Matrix(size_t rows, size_t cols);
    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    ~Matrix();

    size_t getRows() const { return rows_; }
    size_t getCols() const { return cols_; }

    double* operator[](size_t row);
    const double* operator[](size_t row) const;

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;

    Vector operator*(const Vector& vec) const; // Умножение справа
    Vector leftMultiply(const Vector& vec) const; // Умножение слева

    Matrix transpose() const;

    static Matrix generateRandom(size_t rows, size_t cols, double min_val = -10.0, double max_val = 10.0);

    void saveToFile(const std::string& filename) const;

    friend Matrix operator*(double scalar, const Matrix& matrix);
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

private:
    void distributeRows(size_t total_rows, size_t num_threads,
                       std::vector<std::thread>& threads,
                       std::function<void(size_t, size_t)> worker) const;

    // Вспомогательные методы для управления памятью
    void allocateMemory();
    void deallocateMemory();
    void copyData(const Matrix& other);

    size_t rows_;
    size_t cols_;
    double** data_;
};