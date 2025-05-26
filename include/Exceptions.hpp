#pragma once
#include <stdexcept>
#include <exception>

class DimensionMismatchException : public std::exception {
public:
    const char* what() const noexcept override {
        return "Matrix dimensions do not match for the operation.";
    }
};

class InvalidMatrixOperationException : public std::exception {
public:
    const char* what() const noexcept override {
        return "Invalid matrix operation requested.";
    }
};

class OutOfBoundsException : public std::exception {
public:
    const char* what() const noexcept override {
        return "Matrix index out of bounds.";
    }
};
class DivisionByZeroException : public std::exception {
public:
    const char* what() const noexcept override {
        return "Division by zero encountered during matrix operation";
    }
};