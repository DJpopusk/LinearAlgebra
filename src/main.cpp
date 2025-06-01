#include "Matrix.hpp"
#include "Vector.hpp"
#include "ThreadPoolConfig.hpp"
#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

void testDistributeElementsWithDifferentSizes() {
    std::cout << "\n=== distributeElements test===" << std::endl;

    struct TestCase {
        size_t rows, cols, threads;
        std::string description;
    };

    std::vector<TestCase> testCases = {
        {2, 2, 5, "2x2 , 5"},
        {4, 2, 3, "4x2 , 3"},
        {5, 3, 4, "5x3 , 4"},
        {1, 10, 3, "1x10 , 3"},
        {3, 3, 4, "3x3, 4"}
    };

    for (const auto& testCase : testCases) {
        std::cout << "\n" << testCase.description << std::endl;

        ThreadPoolConfig::setNumThreads(testCase.threads);

        size_t total_elements = testCase.rows * testCase.cols;
        size_t effective_threads = std::min(testCase.threads, total_elements);
        size_t base_chunk = total_elements / effective_threads;
        size_t remainder = total_elements % effective_threads;

        std::cout << "Elements: " << total_elements << std::endl;
        std::cout << "Effective threads cnt: " << effective_threads << std::endl;
        std::cout << "Chunk base size: " << base_chunk << std::endl;
        std::cout << "Reminder: " << remainder << std::endl;

        size_t current_element = 0;
        for (size_t t = 0; t < effective_threads; ++t) {
            size_t chunk_size = base_chunk;
            if (remainder > 0 && t >= (effective_threads - remainder)) {
                chunk_size += 1;
            }

            std::cout << "- Thread " << (t + 1) << ": el "
                     << current_element << "-" << (current_element + chunk_size - 1)
                     << " (" << chunk_size << " els)" << std::endl;

            current_element += chunk_size;
        }

        Matrix m1(testCase.rows, testCase.cols);
        Matrix m2(testCase.rows, testCase.cols);

        for (size_t i = 0; i < testCase.rows; ++i) {
            for (size_t j = 0; j < testCase.cols; ++j) {
                m1[i][j] = i * testCase.cols + j + 1;
                m2[i][j] = 1.0;
            }
        }

        Matrix result = m1 + m2;
        std::cout << "✓" << std::endl;
    }
}

Matrix simpleMultiply(const Matrix& a, const Matrix& b) {
    if(a.getCols() != b.getRows()) throw DimensionMismatchException();

    Matrix result(a.getRows(), b.getCols());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < b.getCols(); ++j) {
            double sum = 0.0;
            for(size_t k = 0; k < a.getCols(); ++k) {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

Matrix simpleAdd(const Matrix& a, const Matrix& b) {
    if(a.getRows() != b.getRows() || a.getCols() != b.getCols())
        throw DimensionMismatchException();

    Matrix result(a.getRows(), a.getCols());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < a.getCols(); ++j) {
             result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

Matrix simpleSubtract(const Matrix& a, const Matrix& b) {
    if(a.getRows() != b.getRows() || a.getCols() != b.getCols())
        throw DimensionMismatchException();

    Matrix result(a.getRows(), a.getCols());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < a.getCols(); ++j) {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    return result;
}

Matrix simpleTranspose(const Matrix& a) {
    Matrix result(a.getCols(), a.getRows());
    for(size_t i = 0; i < a.getRows(); ++i) {
        for(size_t j = 0; j < a.getCols(); ++j) {
            result[j][i] = a[i][j];
        }
    }
    return result;
}

bool compareResults(const Matrix& a, const Matrix& b) {
    if(a.getRows() != b.getRows() || a.getCols() != b.getCols()) {
        cout << "Results match: false (different sizes)\n";
        return false;
    }

    bool identical = true;
    for(size_t i = 0; i < a.getRows() && identical; ++i) {
        for(size_t j = 0; j < a.getCols() && identical; ++j) {
            if(std::abs(a[i][j] - b[i][j]) > 1e-6) {
                identical = false;
            }
        }
    }

    cout << "Results match: " << boolalpha << identical << "\n";
    return identical;
}

void testExceptions() {
    cout << "\n=== Testing Exception Handling ===\n";

    // Test DimensionMismatchException
    try {
        Matrix a(2, 3);
        Matrix b(4, 2);
        Matrix result = a + b;
        cout << "ERROR: DimensionMismatchException not thrown!\n";
    } catch (const DimensionMismatchException& e) {
        cout << "✓ Caught DimensionMismatchException: " << e.what() << "\n";
    }

    // Test OutOfBoundsException for Matrix
    try {
        Matrix a(3, 3);
        double val = a[5][2];
        cout << "ERROR: OutOfBoundsException not thrown! val=" << val << "\n";
    } catch (const OutOfBoundsException& e) {
        cout << "✓ Caught OutOfBoundsException (Matrix): " << e.what() << "\n";
    }

    // Test OutOfBoundsException for Vector
    try {
        Vector v(5);
        double val = v[10];
        cout << "ERROR: OutOfBoundsException not thrown! val=" << val << "\n";
    } catch (const OutOfBoundsException& e) {
        cout << "✓ Caught OutOfBoundsException (Vector): " << e.what() << "\n";
    }

    // Test InvalidMatrixOperationException
    try {
        Vector v1(2);
        Vector v2(2);
        Vector cross_result = v1.cross(v2);
        cout << "ERROR: InvalidMatrixOperationException not thrown!\n";
    } catch (const InvalidMatrixOperationException& e) {
        cout << "✓ Caught InvalidMatrixOperationException: " << e.what() << "\n";
    }

    // Test Vector dimension mismatch
    try {
        Vector v1(3);
        Vector v2(5);
        Vector result = v1 + v2;
        cout << "ERROR: DimensionMismatchException not thrown!\n";
    } catch (const DimensionMismatchException& e) {
        cout << "✓ Caught DimensionMismatchException (Vector): " << e.what() << "\n";
    }

    cout << "Exception testing completed!\n";
}

void testPerformance(size_t size) {
    try {
        cout << "\nTesting performance for " << size << "x" << size << " matrices:\n";

        auto start_gen = high_resolution_clock::now();
        Matrix a = Matrix::generateRandom(size, size);
        Matrix b = Matrix::generateRandom(size, size);
        auto end_gen = high_resolution_clock::now();
        cout << "• Generation time: "
             << duration_cast<milliseconds>(end_gen - start_gen).count()
             << " ms\n";

        a.saveToFile("../tests/matrix_a.txt");
        b.saveToFile("../tests/matrix_b.txt");

        auto start_parallel = high_resolution_clock::now();
        Matrix parallel_result = a * b;
        auto end_parallel = high_resolution_clock::now();
        cout << "• Parallel multiply: "
             << duration_cast<milliseconds>(end_parallel - start_parallel).count()
             << " ms\n";
        parallel_result.saveToFile("../tests/cpp_multiply_result.txt");

        auto start_simple = high_resolution_clock::now();
        Matrix simple_result = simpleMultiply(a, b);
        auto end_simple = high_resolution_clock::now();
        cout << "• Simple multiply: "
             << duration_cast<milliseconds>(end_simple - start_simple).count()
             << " ms\n";

        cout << "• ";
        compareResults(parallel_result, simple_result);

        auto start_parallel_add = high_resolution_clock::now();
        Matrix parallel_add = a + b;
        auto end_parallel_add = high_resolution_clock::now();
        cout << "• Parallel add: "
             << duration_cast<milliseconds>(end_parallel_add - start_parallel_add).count()
             << " ms\n";
        parallel_add.saveToFile("../tests/cpp_add_result.txt");

        auto start_simple_add = high_resolution_clock::now();
        Matrix simple_add_result = simpleAdd(a, b);
        auto end_simple_add = high_resolution_clock::now();
        cout << "• Simple add: "
             << duration_cast<milliseconds>(end_simple_add - start_simple_add).count()
             << " ms\n";

        cout << "• ";
        compareResults(parallel_add, simple_add_result);

        auto start_parallel_sub = high_resolution_clock::now();
        Matrix parallel_sub = a - b;
        auto end_parallel_sub = high_resolution_clock::now();
        cout << "• Parallel sub: "
             << duration_cast<milliseconds>(end_parallel_sub - start_parallel_sub).count()
             << " ms\n";
        parallel_sub.saveToFile("../tests/cpp_sub_result.txt");

        auto start_simple_sub = high_resolution_clock::now();
        Matrix simple_sub_result = simpleSubtract(a, b);
        auto end_simple_sub = high_resolution_clock::now();
        cout << "• Simple sub: "
             << duration_cast<milliseconds>(end_simple_sub - start_simple_sub).count()
             << " ms\n";

        cout << "• ";
        compareResults(parallel_sub, simple_sub_result);

        auto start_parallel_transpose = high_resolution_clock::now();
        Matrix parallel_transpose = a.transpose();
        auto end_parallel_transpose = high_resolution_clock::now();
        cout << "• Parallel transpose: "
             << duration_cast<milliseconds>(end_parallel_transpose - start_parallel_transpose).count()
             << " ms\n";
        parallel_transpose.saveToFile("../tests/cpp_transpose_result.txt");

        auto start_simple_transpose = high_resolution_clock::now();
        Matrix simple_transpose_result = simpleTranspose(a);
        auto end_simple_transpose = high_resolution_clock::now();
        cout << "• Simple transpose: "
             << duration_cast<milliseconds>(end_simple_transpose - start_simple_transpose).count()
             << " ms\n";

        cout << "• ";
        compareResults(parallel_transpose, simple_transpose_result);

    } catch(const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
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
        testExceptions();
        testPerformance(500);
        testPerformance(1000);

        cout << "\n=== Files saved for Python verification ===\n";
        cout << "Run 'python3 verify_results.py' to verify results with NumPy\n";

    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    try {
        testDistributeElementsWithDifferentSizes();

        std::cout << "\n=== Все тесты завершены ===" << std::endl;

    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}