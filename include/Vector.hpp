#pragma once
#include <cstddef>
#include <iostream>
#include <thread>
#include <vector>
#include <functional>
#include "Exceptions.hpp"
#include "ThreadPoolConfig.hpp"

class Vector {
public:
    Vector(size_t size);
    Vector(const Vector& other);
    Vector& operator=(const Vector& other);
    ~Vector();

    size_t size() const { return size_; }

    double& operator[](size_t idx);
    const double& operator[](size_t idx) const;

    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;

    double dot(const Vector& other) const;
    Vector cross(const Vector& other) const; // Только для размерности 3
    Vector transpose() const;

    static Vector generateRandom(size_t size, double min_val = -10.0, double max_val = 10.0);

    void saveToFile(const std::string& filename) const;

    friend Vector operator*(double scalar, const Vector& vector);
    friend std::ostream& operator<<(std::ostream& os, const Vector& vector);

private:
    void distributeElements(size_t total_elements, size_t num_threads,
                           std::vector<std::thread>& threads,
                           std::function<void(size_t, size_t)> worker) const;
    
    size_t size_;
    double* data_;
};