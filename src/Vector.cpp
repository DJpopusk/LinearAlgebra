#include "Vector.hpp"
#include <cstring>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <random>
#include <fstream>

Vector::Vector(size_t size) : size_(size), data_(new double[size]()) {}

Vector::Vector(const Vector& other) : size_(other.size_), data_(new double[other.size_]) {
    std::memcpy(data_, other.data_, size_ * sizeof(double));
}

Vector& Vector::operator=(const Vector& other) {
    if (this != &other) {
        delete[] data_;
        size_ = other.size_;
        data_ = new double[size_];
        std::memcpy(data_, other.data_, size_ * sizeof(double));
    }
    return *this;
}

Vector::~Vector() {
    delete[] data_;
}

double& Vector::operator[](size_t idx) {
    if (idx >= size_) throw OutOfBoundsException();
    return data_[idx];
}

const double& Vector::operator[](size_t idx) const {
    if (idx >= size_) throw OutOfBoundsException();
    return data_[idx];
}

void Vector::distributeElements(size_t total_elements, size_t num_threads,
                               std::vector<std::thread>& threads,
                               std::function<void(size_t, size_t)> worker) const {
    const size_t chunk = (total_elements + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end = std::min(start + chunk, total_elements);
        if (start < end) {
            threads.emplace_back(worker, start, end);
        }
    }
    for (std::thread& t : threads) t.join();
    threads.clear();
}

Vector Vector::operator+(const Vector& other) const {
    if (size_ != other.size_) throw DimensionMismatchException();

    Vector result(size_);
    const size_t num_threads = std::min(size_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
    };

    distributeElements(size_, num_threads, threads, worker);
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    if (size_ != other.size_) throw DimensionMismatchException();

    Vector result(size_);
    const size_t num_threads = std::min(size_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
    };

    distributeElements(size_, num_threads, threads, worker);
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(size_);
    const size_t num_threads = std::min(size_, ThreadPoolConfig::getNumThreads());
    std::vector<std::thread> threads;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            result.data_[i] = data_[i] * scalar;
        }
    };

    distributeElements(size_, num_threads, threads, worker);
    return result;
}

double Vector::dot(const Vector& other) const {
    if (size_ != other.size_) throw DimensionMismatchException();

    const size_t num_threads = std::min(size_, ThreadPoolConfig::getNumThreads());
    std::vector<double> partial_sums(num_threads, 0.0);
    std::vector<std::thread> threads;

    const size_t chunk = (size_ + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end = std::min(start + chunk, size_);
        if (start < end) {
            threads.emplace_back([&, t, start, end]() {
                double local_sum = 0.0;
                for (size_t i = start; i < end; ++i) {
                    local_sum += data_[i] * other.data_[i];
                }
                partial_sums[t] = local_sum;
            });
        }
    }

    for (std::thread& t : threads) t.join();

    double total_sum = 0.0;
    for (double sum : partial_sums) {
        total_sum += sum;
    }

    return total_sum;
}

Vector Vector::cross(const Vector& other) const {
    if (size_ != 3 || other.size_ != 3) {
        throw InvalidMatrixOperationException();
    }

    Vector result(3);
    result.data_[0] = data_[1] * other.data_[2] - data_[2] * other.data_[1];
    result.data_[1] = data_[2] * other.data_[0] - data_[0] * other.data_[2];
    result.data_[2] = data_[0] * other.data_[1] - data_[1] * other.data_[0];

    return result;
}

Vector Vector::transpose() const {
    return Vector(*this);
}

Vector Vector::generateRandom(size_t size, double min_val, double max_val) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_val, max_val);

    Vector v(size);
    for(size_t i = 0; i < size; ++i) {
        v.data_[i] = dis(gen);
    }
    return v;
}

void Vector::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    file << size_ << std::endl;
    for (size_t i = 0; i < size_; ++i) {
        file << std::fixed << std::setprecision(6) << data_[i];
        if (i < size_ - 1) file << " ";
    }
    file << std::endl;
}

Vector operator*(double scalar, const Vector& vector) {
    return vector * scalar;
}

std::ostream& operator<<(std::ostream& os, const Vector& vector) {
    os << "[";
    for (size_t i = 0; i < vector.size_; ++i) {
        os << std::setw(8) << std::fixed << std::setprecision(2) << vector.data_[i];
        if (i < vector.size_ - 1) os << ", ";
    }
    os << "]";
    return os;
}