#pragma once
#include <fstream>
#include <string>
#include "Matrix.hpp"
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
};