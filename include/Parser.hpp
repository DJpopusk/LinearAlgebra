#pragma once
#include <fstream>
#include <string>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Exceptions.hpp"

class MatrixParser {
public:
    static Matrix readMatrixFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw ParserException("Failed to open file: " + filename);
        }

        size_t rows = 0, cols = 0;
        file >> rows >> cols;
        if (!file.good() || rows == 0 || cols == 0) {
            throw ParserException("Invalid matrix dimensions in file: " + filename);
        }

        Matrix mat(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                double value = 0.0;
                if (!(file >> value)) {
                    throw ParserException(
                        "Failed to read matrix element at row "
                        + std::to_string(i)
                        + ", col "
                        + std::to_string(j)
                        + " in file: " + filename
                    );
                }
                mat(i, j) = value;
            }
        }
        return mat;
    }

    static Vector readVectorFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw ParserException("Failed to open file: " + filename);
        }

        size_t size = 0;
        file >> size;
        if (!file.good() || size == 0) {
            throw ParserException("Invalid vector size in file: " + filename);
        }

        Vector vec(size);
        for (size_t i = 0; i < size; ++i) {
            double value = 0.0;
            if (!(file >> value)) {
                throw ParserException(
                    "Failed to read vector element at index "
                    + std::to_string(i)
                    + " in file: " + filename
                );
            }
            vec(i) = value;
        }
        return vec;
    }
};