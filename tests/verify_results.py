import numpy as np
import os

def load_matrix_from_file(filename):
    if not os.path.exists(filename):
        print(f"File {filename} not found!")
        return None

    with open(filename, 'r') as f:
        lines = f.readlines()

    rows, cols = map(int, lines[0].strip().split())

    matrix = np.zeros((rows, cols))
    for i in range(rows):
        values = list(map(float, lines[i+1].strip().split()))
        matrix[i] = values

    return matrix

def verify_operations():
    print("=== Python NumPy Verification ===")

    A = load_matrix_from_file("matrix_a.txt")
    B = load_matrix_from_file("matrix_b.txt")

    if A is None or B is None:
        print("Failed to load input matrices!")
        return

    print(f"Loaded matrices: A{A.shape}, B{B.shape}")

    cpp_multiply = load_matrix_from_file("cpp_multiply_result.txt")
    if cpp_multiply is not None:
        numpy_multiply = A @ B
        diff = np.abs(cpp_multiply - numpy_multiply)
        max_diff = np.max(diff)
        print(f"✓ Matrix multiplication - Max difference: {max_diff:.2e}")
        if max_diff < 1e-3:
            print("  PASSED: Results match within tolerance")
        else:
            print("  FAILED: Results differ significantly")

    cpp_add = load_matrix_from_file("cpp_add_result.txt")
    if cpp_add is not None:
        numpy_add = A + B
        diff = np.abs(cpp_add - numpy_add)
        max_diff = np.max(diff)
        print(f"✓ Matrix addition - Max difference: {max_diff:.2e}")
        if max_diff < 2e-6:
            print("  PASSED: Results match within tolerance")
        else:
            print("  FAILED: Results differ significantly")

    cpp_sub = load_matrix_from_file("cpp_sub_result.txt")
    if cpp_sub is not None:
        numpy_sub = A - B
        diff = np.abs(cpp_sub - numpy_sub)
        max_diff = np.max(diff)
        print(f"✓ Matrix subtraction - Max difference: {max_diff:.2e}")
        if max_diff < 2e-6:
            print("  PASSED: Results match within tolerance")
        else:
            print("  FAILED: Results differ significantly")

    cpp_transpose = load_matrix_from_file("cpp_transpose_result.txt")
    if cpp_transpose is not None:
        numpy_transpose = A.T
        diff = np.abs(cpp_transpose - numpy_transpose)
        max_diff = np.max(diff)
        print(f"✓ Matrix transpose - Max difference: {max_diff:.2e}")
        if max_diff < 1e-10:
            print("  PASSED: Results match within tolerance")
        else:
            print("  FAILED: Results differ significantly")

    print("\nVerification completed!")

if __name__ == "__main__":
    verify_operations()