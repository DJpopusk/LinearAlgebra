#include "Matrix.hpp"
#include <thread>
#include <vector>


Matrix::Matrix(size_t rows, size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols) {}

double& Matrix::operator()(size_t row, size_t col) {
    if (row >= rows_ || col >= cols_) throw OutOfBoundsException();
    return data_[row * cols_ + col];
}

const double& Matrix::operator()(size_t row, size_t col) const {
    if (row >= rows_ || col >= cols_) throw OutOfBoundsException();
    return data_[row * cols_ + col];
}

void Matrix::distribute_rows(size_t total_rows, size_t num_threads,
                            std::vector<std::thread>& threads,
                            std::function<void(size_t, size_t)> worker) const {
    const size_t chunk = (total_rows + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end = std::min(start + chunk, total_rows);
        if (start < end) {
            threads.emplace_back(worker, start, end);
        }
    }
    for (auto& t : threads) t.join();
    threads.clear();
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw DimensionMismatchException();
    Matrix result(rows_, cols_);
    const size_t num_threads = std::min(rows_, static_cast<size_t>(std::thread::hardware_concurrency()));
    std::vector<std::thread> threads;
    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
    };
    distribute_rows(rows_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw DimensionMismatchException();
    Matrix result(rows_, cols_);
    const size_t num_threads = std::min(rows_, static_cast<size_t>(std::thread::hardware_concurrency()));
    std::vector<std::thread> threads;
    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
    };
    distribute_rows(rows_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) throw DimensionMismatchException();
    Matrix result(rows_, other.cols_);
    const size_t num_threads = std::min(rows_, static_cast<size_t>(std::thread::hardware_concurrency()));
    std::vector<std::thread> threads;
    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < cols_; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
    };
    distribute_rows(rows_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(rows_, cols_);
    const size_t num_threads = std::min(rows_, static_cast<size_t>(std::thread::hardware_concurrency()));
    std::vector<std::thread> threads;
    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }
    };
    distribute_rows(rows_, num_threads, threads, worker);
    return result;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    return matrix * scalar;
}
