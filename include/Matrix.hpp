#pragma once
#include <cstddef>
#include <iostream>
#include <thread>
#include <vector>
#include <functional>
#include "Exceptions.hpp"
#include "ThreadPoolConfig.hpp"

class Vector; // Forward declaration

class Matrix {
public:
    Matrix(size_t rows, size_t cols);
    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    ~Matrix();

    size_t getRows() const { return rows_; }
    size_t getCols() const { return cols_; }

    double& operator()(size_t row, size_t col);
    const double& operator()(size_t row, size_t col) const;

    // Базовые арифметические операции
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;

    // Операции с векторами
    Vector operator*(const Vector& vec) const; // Умножение справа
    Vector leftMultiply(const Vector& vec) const; // Умножение слева

    // Специальные операции
    Matrix transpose() const;

    // Операции для повышенного коэффициента
    double determinant() const;
    Matrix inverse() const;
    double maxEigenvalue() const;
    Vector solve(const Vector& b) const; // Решение СЛАУ

    // Операции для максимального коэффициента
    std::vector<double> eigenvalues() const;

    // Дружественные функции
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

    Matrix luDecomposition(std::vector<size_t>& permutation) const;
    double determinantLU(const Matrix& lu, const std::vector<size_t>& perm) const;
    Vector forwardSubstitution(const Matrix& L, const Vector& b) const;
    Vector backwardSubstitution(const Matrix& U, const Vector& y) const;

    size_t rows_;
    size_t cols_;
    double** data_; // Двумерный массив указателей
};