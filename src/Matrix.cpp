#include "Matrix.hpp"
#include "Vector.hpp"
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>

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
        data_[i] = new double[cols_]();
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
    const size_t total_elements = rows_ * cols_;
    const size_t num_threads = std::min(total_elements, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem / cols_;
            size_t j = elem % cols_;
            data_[i][j] = other.data_[i][j];
        }
    };

    distributeElements(total_elements, num_threads, threads, worker);
}

double* Matrix::operator[](size_t row) {
    if (row >= rows_) throw OutOfBoundsException();
    return data_[row];
}

const double* Matrix::operator[](size_t row) const {
    if (row >= rows_) throw OutOfBoundsException();
    return data_[row];
}

void Matrix::distributeElements(size_t total_elements, size_t num_threads,
                               std::vector<std::thread>& threads,
                               std::function<void(size_t, size_t)> worker) const {
    if (total_elements == 0) {
        return; // Нет работы для выполнения
    }

    // Ограничиваем количество потоков количеством элементов
    const size_t effective_threads = std::min(num_threads, total_elements);

    if (effective_threads == 0) {
        return;
    }

    // Базовое количество элементов на поток (округление вниз)
    const size_t base_chunk = total_elements / effective_threads;
    // Количество оставшихся элементов для распределения
    const size_t remainder = total_elements % effective_threads;

    size_t current_element = 0;

    for (size_t t = 0; t < effective_threads; ++t) {
        size_t chunk_size = base_chunk;

        // Если есть оставшиеся элементы, добавляем их к последним потокам
        if (remainder > 0 && t >= (effective_threads - remainder)) {
            chunk_size += 1;
        }

        const size_t start_element = current_element;
        const size_t end_element = start_element + chunk_size;

        if (start_element < total_elements && chunk_size > 0) {
            threads.emplace_back(worker, start_element, end_element);
        }

        current_element += chunk_size;
    }

    // Ожидаем завершения всех потоков
    for (std::thread& t : threads) {
        t.join();
    }
    threads.clear();
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw DimensionMismatchException();

    Matrix result(rows_, cols_);
    const size_t total_elements = rows_ * cols_;
    const size_t num_threads = std::min(total_elements, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem / cols_;
            size_t j = elem % cols_;
            result.data_[i][j] = data_[i][j] + other.data_[i][j];
        }
    };

    distributeElements(total_elements, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw DimensionMismatchException();

    Matrix result(rows_, cols_);
    const size_t total_elements = rows_ * cols_;
    const size_t num_threads = std::min(total_elements, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem / cols_;
            size_t j = elem % cols_;
            result.data_[i][j] = data_[i][j] - other.data_[i][j];
        }
    };

    distributeElements(total_elements, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) throw DimensionMismatchException();

    Matrix result(rows_, other.cols_);
    const size_t total_elements = rows_ * other.cols_;
    const size_t num_threads = std::min(total_elements, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem / other.cols_;
            size_t j = elem % other.cols_;
            double sum = 0.0;
            for (size_t k = 0; k < cols_; ++k) {
                sum += data_[i][k] * other.data_[k][j];
            }
            result.data_[i][j] = sum;
        }
    };

    distributeElements(total_elements, num_threads, threads, worker);
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(rows_, cols_);
    const size_t total_elements = rows_ * cols_;
    const size_t num_threads = std::min(total_elements, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem / cols_;
            size_t j = elem % cols_;
            result.data_[i][j] = data_[i][j] * scalar;
        }
    };

    distributeElements(total_elements, num_threads, threads, worker);
    return result;
}

Vector Matrix::operator*(const Vector& vec) const {
    if (cols_ != vec.size()) throw DimensionMismatchException();

    Vector result(rows_);
    const size_t num_threads = std::min(rows_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem; // для Vector умножения элемент = строка
            double sum = 0.0;
            for (size_t j = 0; j < cols_; ++j) {
                sum += data_[i][j] * vec[j];
            }
            result[i] = sum;
        }
    };

    distributeElements(rows_, num_threads, threads, worker);
    return result;
}

Vector Matrix::leftMultiply(const Vector& vec) const {
    if (rows_ != vec.size()) throw DimensionMismatchException();

    Vector result(cols_);
    const size_t num_threads = std::min(cols_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t j = elem; // для Vector умножения элемент = столбец
            double sum = 0.0;
            for (size_t i = 0; i < rows_; ++i) {
                sum += vec[i] * data_[i][j];
            }
            result[j] = sum;
        }
    };

    distributeElements(cols_, num_threads, threads, worker);
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    const size_t total_elements = rows_ * cols_;
    const size_t num_threads = std::min(total_elements, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start_elem, size_t end_elem) {
        for (size_t elem = start_elem; elem < end_elem; ++elem) {
            size_t i = elem / cols_;
            size_t j = elem % cols_;
            result.data_[j][i] = data_[i][j];
        }
    };

    distributeElements(total_elements, num_threads, threads, worker);
    return result;
}

Matrix Matrix::generateRandom(size_t rows, size_t cols, double min_val, double max_val) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_val, max_val);

    Matrix m(rows, cols);
    for(size_t i = 0; i < rows; ++i) {
        for(size_t j = 0; j < cols; ++j) {
            m.data_[i][j] = dis(gen);
        }
    }
    return m;
}

void Matrix::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    file << rows_ << " " << cols_ << std::endl;
    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < cols_; ++j) {
            file << std::fixed << std::setprecision(6) << data_[i][j];
            if (j < cols_ - 1) file << " ";
        }
        file << std::endl;
    }
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