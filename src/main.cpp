#include "Matrix.hpp"
#include "Parser.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;


Matrix generate_random_matrix(size_t rows, size_t cols) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-10.0, 10.0);

    Matrix m(rows, cols);
    for(size_t i = 0; i < rows; ++i) {
        for(size_t j = 0; j < cols; ++j) {
            m(i, j) = dis(gen);
        }
    }
    return m;
}

Matrix simple_multiply(const Matrix& a, const Matrix& b) {
    if(a.getCols() != b.getRows()) throw DimensionMismatchException();

    Matrix result(a.getRows(), b.getCols());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < b.getCols(); ++j) {
            double sum = 0.0;
            for(size_t k = 0; k < a.getCols(); ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

Matrix simple_add(const Matrix& a, const Matrix& b) {
    if(a.getRows() != b.getRows() || a.getCols() != b.getCols()) throw DimensionMismatchException();

    Matrix result(a.getRows(), a.getCols());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < a.getCols(); ++j) {
             result(i,j) = a(i, j) + b(i, j);
        }
    }
    return result;
}

Matrix simple_subtract(const Matrix& a, const Matrix& b) {
    if(a.getRows() != b.getRows() || a.getCols() != b.getCols()) throw DimensionMismatchException();

    Matrix result(a.getRows(), a.getCols());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < a.getCols(); ++j) {
            result(i,j) = a(i, j) - b(i, j);
        }
    }
    return result;
}

bool same_res(const Matrix& a, const Matrix& b) {
    if(a.getRows() != b.getRows() || a.getCols() != b.getCols()) {
        cout << "Results match: false (different sizes)\n";
        return false;
    }

    bool identical = true;
    for(size_t i = 0; i < a.getRows() && identical; ++i) {
        for(size_t j = 0; j < a.getCols() && identical; ++j) {
            if(std::abs(a(i, j) - b(i, j)) > 1e-6) {
                identical = false;
            }
        }
    }

    cout << "Results match: " << boolalpha << identical << "\n";
    return identical;
}

void test_performance(size_t size) {
    try {
        cout << "\nTesting performance for " << size << "x" << size << " matrices:\n";

        auto start_gen = high_resolution_clock::now();
        Matrix a = generate_random_matrix(size, size);
        Matrix b = generate_random_matrix(size, size);
        auto end_gen = high_resolution_clock::now();
        cout << "Generation time: "
             << duration_cast<milliseconds>(end_gen - start_gen).count()
             << " ms\n";

        auto start_parallel = high_resolution_clock::now();
        Matrix parallel_result = a * b;
        auto end_parallel = high_resolution_clock::now();
        cout << "Parallel multiply: "
             << duration_cast<milliseconds>(end_parallel - start_parallel).count()
             << " ms\n";

        auto start_simple = high_resolution_clock::now();
        Matrix simple_result = simple_multiply(a, b);
        auto end_simple = high_resolution_clock::now();
        cout << "Simple multiply: "
             << duration_cast<milliseconds>(end_simple - start_simple).count()
             << " ms\n";

        same_res(parallel_result, simple_result);

        auto start_parallel_add = high_resolution_clock::now();
        Matrix parallel_add = a + b;
        auto end_parallel_add = high_resolution_clock::now();
        cout << "Parallel add: "
             << duration_cast<milliseconds>(end_parallel_add - start_parallel_add).count()
             << " ms\n";

        auto start_simple_add = high_resolution_clock::now();
        Matrix simple_add_result = simple_add(a, b);
        auto end_simple_add = high_resolution_clock::now();
        cout << "Simple add: "
             << duration_cast<milliseconds>(end_simple_add - start_simple_add).count()
             << " ms\n";

        same_res(parallel_add, simple_add_result);

        auto start_parallel_sub = high_resolution_clock::now();
        Matrix parallel_sub = a - b;
        auto end_parallel_sub = high_resolution_clock::now();
        cout << "Parallel sub: "
             << duration_cast<milliseconds>(end_parallel_sub - start_parallel_sub).count()
             << " ms\n";

        auto start_simple_sub = high_resolution_clock::now();
        Matrix simple_sub_result = simple_subtract(a, b);
        auto end_simple_sub = high_resolution_clock::now();
        cout << "Simple sub: "
             << duration_cast<milliseconds>(end_simple_sub - start_simple_sub).count()
             << " ms\n";

        same_res(parallel_sub, simple_sub_result);

    } catch(const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
}

int main() {
    try {

        Matrix A = MatrixParser::readMatrixFromFile("/Users/phonkyponky/CLionProjects/LinearAlgebra/tests/matrixA.txt");

        Matrix B = MatrixParser::readMatrixFromFile("/Users/phonkyponky/CLionProjects/LinearAlgebra/tests/matrixB.txt");

        Matrix sumAB = A + B;
        Matrix diffAB = A - B;
        Matrix mulAB = A * B;
        Matrix scalarMul = A * 2.5;

        std::cout << "Matrix A:\n" << A << std::endl;
        std::cout << "Matrix B:\n" << B << std::endl;

        std::cout << "A + B:\n" << sumAB << std::endl;
        std::cout << "A - B:\n" << diffAB << std::endl;
        std::cout << "A * B:\n" << mulAB << std::endl;
        std::cout << "A * 2.5:\n" << scalarMul << std::endl;
    }
    catch (const std::exception& ex) {
        std::cerr << "Ошибка: " << ex.what() << std::endl;
    }
    // Correct operations

    return 0;
}