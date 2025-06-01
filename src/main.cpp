#include "Matrix.hpp"
#include "Vector.hpp"
#include "Parser.hpp"
#include "ThreadPoolConfig.hpp"
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

void testMatrixOperations() {
    cout << "\n=== Testing Matrix Operations with double** ===\n";

    Matrix A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;

    Matrix B(2, 3);
    B(0, 0) = 7; B(0, 1) = 8; B(0, 2) = 9;
    B(1, 0) = 10; B(1, 1) = 11; B(1, 2) = 12;

    cout << "Matrix A:\n" << A << endl;
    cout << "Matrix B:\n" << B << endl;

    Matrix sum = A + B;
    cout << "A + B:\n" << sum << endl;

    Matrix diff = A - B;
    cout << "A - B:\n" << diff << endl;

    Matrix scaled = A * 2.5;
    cout << "A * 2.5:\n" << scaled << endl;

    Matrix AT = A.transpose();
    cout << "A^T:\n" << AT << endl;
}

int main(int argc, char* argv[]) {
    size_t num_threads = thread::hardware_concurrency();

    if (argc > 1) {
        num_threads = stoul(argv[1]);
        cout << "Using " << num_threads << " threads (from command line)\n";
    } else {
        cout << "Enter number of threads (default " << num_threads << "): ";
        string input;
        getline(cin, input);
        if (!input.empty()) {
            num_threads = stoul(input);
        }
        cout << "Using " << num_threads << " threads\n";
    }

    ThreadPoolConfig::setNumThreads(num_threads);

    try {
        testMatrixOperations();

        cout << "\n=== Performance Testing with double** ===\n";

        Matrix A = generate_random_matrix(100, 100);
        Matrix B = generate_random_matrix(100, 100);

        auto start = high_resolution_clock::now();
        Matrix C = A * B;
        auto end = high_resolution_clock::now();

        cout << "Matrix multiplication (100x100): "
             << duration_cast<milliseconds>(end - start).count()
             << " ms\n";

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}