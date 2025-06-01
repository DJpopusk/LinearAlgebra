#include "Matrix.hpp"
#include "Vector.hpp"
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <random>

Matrix::Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
    allocateMemory();
}

Matrix::Matrix(const Matrix& other) : rows_(other.rows_), cols_(other.cols_) {
    allocateMemory();
    copyData(other);
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        deallocateMemory();
        rows_ = other.rows_;
        cols_ = other.cols_;
        allocateMemory();
        copyData(other);
    }
    return *this;
}

Matrix::~Matrix() {
    deallocateMemory();
}

void Matrix::allocateMemory() {
    data_ = new double*[rows_];
    for (size_t i = 0; i < rows_; ++i) {
        data_[i] = new double[cols_]();  // Инициализация нулями
    }
}

void Matrix::deallocateMemory() {
    if (data_) {
        for (size_t i = 0; i < rows_; ++i) {
            delete[] data_[i];
        }
        delete[] data_;
        data_ = nullptr;
    }
}

void Matrix::copyData(const Matrix& other) {
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                data_[i][j] = other.data_[i][j];
            }
        }
    };

    distributeRows(rows_, num_threads, threads, worker);
}

double& Matrix::operator()(size_t row, size_t col) {
    if (row >= rows_ || col >= cols_) throw OutOfBoundsException();
    return data_[row][col];
}

const double& Matrix::operator()(size_t row, size_t col) const {
    if (row >= rows_ || col >= cols_) throw OutOfBoundsException();
    return data_[row][col];
}

void Matrix::distributeRows(size_t total_rows, size_t num_threads,
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
    for (std::thread& t : threads) t.join();
    threads.clear();
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw DimensionMismatchException();

    Matrix result(rows_, cols_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data_[i][j] = data_[i][j] + other.data_[i][j];
            }
        }
    };

    distributeRows(rows_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw DimensionMismatchException();

    Matrix result(rows_, cols_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data_[i][j] = data_[i][j] - other.data_[i][j];
            }
        }
    };

    distributeRows(rows_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) throw DimensionMismatchException();

    Matrix result(rows_, other.cols_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < other.cols_; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < cols_; ++k) {
                    sum += data_[i][k] * other.data_[k][j];
                }
                result.data_[i][j] = sum;
            }
        }
    };

    distributeRows(rows_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(rows_, cols_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data_[i][j] = data_[i][j] * scalar;
            }
        }
    };

    distributeRows(rows_, num_threads, threads, worker);
    return result;
}

Vector Matrix::operator*(const Vector& vec) const {
    if (cols_ != vec.size()) throw DimensionMismatchException();

    Vector result(rows_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    const size_t chunk = (rows_ + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end = std::min(start + chunk, rows_);
        if (start < end) {
            threads.emplace_back([&, start, end]() {
                for (size_t i = start; i < end; ++i) {
                    double sum = 0.0;
                    for (size_t j = 0; j < cols_; ++j) {
                        sum += data_[i][j] * vec(j);
                    }
                    result(i) = sum;
                }
            });
        }
    }

    for (std::thread& t : threads) t.join();
    return result;
}

Vector Matrix::leftMultiply(const Vector& vec) const {
    if (rows_ != vec.size()) throw DimensionMismatchException();

    Vector result(cols_);
    const size_t num_threads = std::min(cols_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    const size_t chunk = (cols_ + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end = std::min(start + chunk, cols_);
        if (start < end) {
            threads.emplace_back([&, start, end]() {
                for (size_t j = start; j < end; ++j) {
                    double sum = 0.0;
                    for (size_t i = 0; i < rows_; ++i) {
                        sum += vec(i) * data_[i][j];
                    }
                    result(j) = sum;
                }
            });
        }
    }

    for (std::thread& t : threads) t.join();
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result.data_[j][i] = data_[i][j];
            }
        }
    };

    distributeRows(rows_, num_threads, threads, worker);
    return result;
}

// Простая реализация определителя через LU-разложение
double Matrix::determinant() const {
    if (rows_ != cols_) throw InvalidMatrixOperationException();

    std::vector<size_t> permutation;
    Matrix lu = luDecomposition(permutation);
    return determinantLU(lu, permutation);
}

// Упрощенная реализация LU-разложения с double**
Matrix Matrix::luDecomposition(std::vector<size_t>& permutation) const {
    if (rows_ != cols_) throw InvalidMatrixOperationException();

    Matrix result(*this);
    permutation.resize(rows_);
    for (size_t i = 0; i < rows_; ++i) {
        permutation[i] = i;
    }

    for (size_t k = 0; k < rows_ - 1; ++k) {
        // Найти главный элемент
        size_t max_row = k;
        for (size_t i = k + 1; i < rows_; ++i) {
            if (std::abs(result.data_[i][k]) > std::abs(result.data_[max_row][k])) {
                max_row = i;
            }
        }

        // Обменять строки
        if (max_row != k) {
            std::swap(result.data_[k], result.data_[max_row]);
            std::swap(permutation[k], permutation[max_row]);
        }

        if (std::abs(result.data_[k][k]) < 1e-12) {
            throw DivisionByZeroException();
        }

        // Исключение Гаусса
        for (size_t i = k + 1; i < rows_; ++i) {
            result.data_[i][k] /= result.data_[k][k];
            for (size_t j = k + 1; j < cols_; ++j) {
                result.data_[i][j] -= result.data_[i][k] * result.data_[k][j];
            }
        }
    }

    return result;
}

double Matrix::determinantLU(const Matrix& lu, const std::vector<size_t>& perm) const {
    double det = 1.0;
    size_t swaps = 0;

    for (size_t i = 0; i < rows_; ++i) {
        det *= lu.data_[i][i];
        if (perm[i] != i) swaps++;
    }

    return (swaps % 2 == 0) ? det : -det;
}

// Упрощенная реализация обратной матрицы
Matrix Matrix::inverse() const {
    if (rows_ != cols_) throw InvalidMatrixOperationException();

    std::vector<size_t> permutation;
    Matrix lu = luDecomposition(permutation);
    Matrix result(rows_, cols_);

    // Решить систему для каждого столбца единичной матрицы
    for (size_t j = 0; j < cols_; ++j) {
        Vector e(rows_);
        e(j) = 1.0;

        // Применить перестановки
        Vector b(rows_);
        for (size_t i = 0; i < rows_; ++i) {
            b(i) = e(permutation[i]);
        }

        Vector y = forwardSubstitution(lu, b);
        Vector x = backwardSubstitution(lu, y);

        for (size_t i = 0; i < rows_; ++i) {
            result.data_[i][j] = x(i);
        }
    }

    return result;
}

Vector Matrix::forwardSubstitution(const Matrix& L, const Vector& b) const {
    Vector x(b.size());

    for (size_t i = 0; i < rows_; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < i; ++j) {
            sum += L.data_[i][j] * x(j);
        }
        x(i) = b(i) - sum;
    }

    return x;
}

Vector Matrix::backwardSubstitution(const Matrix& U, const Vector& y) const {
    Vector x(y.size());

    for (int i = static_cast<int>(rows_) - 1; i >= 0; --i) {
        double sum = 0.0;
        for (size_t j = i + 1; j < cols_; ++j) {
            sum += U.data_[i][j] * x(j);
        }
        x(i) = (y(i) - sum) / U.data_[i][i];
    }

    return x;
}

Vector Matrix::solve(const Vector& b) const {
    if (rows_ != b.size()) throw DimensionMismatchException();

    std::vector<size_t> permutation;
    Matrix lu = luDecomposition(permutation);

    // Применить перестановки к правой части
    Vector pb(rows_);
    for (size_t i = 0; i < rows_; ++i) {
        pb(i) = b(permutation[i]);
    }

    Vector y = forwardSubstitution(lu, pb);
    return backwardSubstitution(lu, y);
}

// Простая реализация степенного метода для максимального собственного значения
double Matrix::maxEigenvalue() const {
    if (rows_ != cols_) throw InvalidMatrixOperationException();

    const double tolerance = 1e-8;
    const size_t max_iterations = 1000;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    Vector x(rows_);
    for (size_t i = 0; i < rows_; ++i) {
        x(i) = dis(gen);
    }

    double eigenvalue = 0.0;
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        Vector y = (*this) * x;

        // Нормализация
        double norm = 0.0;
        for (size_t i = 0; i < rows_; ++i) {
            norm += y(i) * y(i);
        }
        norm = std::sqrt(norm);

        if (norm < 1e-12) throw InvalidMatrixOperationException();

        for (size_t i = 0; i < rows_; ++i) {
            y(i) /= norm;
        }

        double new_eigenvalue = x.dot((*this) * x) / x.dot(x);

        if (std::abs(new_eigenvalue - eigenvalue) < tolerance) {
            return new_eigenvalue;
        }

        eigenvalue = new_eigenvalue;
        x = y;
    }

    return eigenvalue;
}

// Заглушка для всех собственных значений (сложная операция)
std::vector<double> Matrix::eigenvalues() const {
    if (rows_ != cols_) throw InvalidMatrixOperationException();

    // Упрощенная реализация - возвращаем только максимальное собственное значение
    // В полной реализации здесь должен быть QR-алгоритм или метод Якоби
    std::vector<double> result;
    result.push_back(maxEigenvalue());

    return result;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    return matrix * scalar;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (size_t i = 0; i < matrix.rows_; i++) {
        for (size_t j = 0; j < matrix.cols_; j++) {
            os << std::setw(8) << std::fixed << std::setprecision(2)
               << matrix.data_[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}