#pragma once
#include <fstream>
#include <stdexcept>
#include <string>
#include "Matrix.hpp"

class MatrixParser {
public:
    // Функция чтения матрицы из текстового файла.
    // Ожидается, что первая строка файла содержит количество строк (rows),
    // вторая строка — количество столбцов (cols),
    // далее построчно располагаются элементы матрицы.
    static Matrix readMatrixFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        size_t rows = 0, cols = 0;
        file >> rows >> cols;
        if (!file.good() || rows == 0 || cols == 0) {
            throw std::runtime_error("Invalid matrix dimensions in file: " + filename);
        }

        Matrix mat(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                double value = 0.0;
                file >> value;
                if (!file.good()) {
                    throw std::runtime_error("Failed to read matrix element at row " 
                                             + std::to_string(i) 
                                             + ", col " 
                                             + std::to_string(j));
                }
                mat(i, j) = value;
            }
        }
        return mat;
    }
};