#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>

#define MAX_SIZE 10

// Function to multiply matrix by vector
void matrix_vector_mult(double A[][MAX_SIZE], int rows, int cols, double v[], double result[]) {
    for (int i = 0; i < rows; i++) {
        result[i] = 0;
        for (int j = 0; j < cols; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
}

// Function to transpose a real matrix
void transpose(double B[][MAX_SIZE], int rows, int cols, double result[][MAX_SIZE]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = B[i][j];
        }
    }
}

// Function to compute conjugate transpose of complex matrix
void conjugate_transpose(double complex B[][MAX_SIZE], int rows, int cols, double complex result[][MAX_SIZE]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = conj(B[i][j]);
        }
    }
}

// Function to compute unconjugated transpose of complex matrix
void unconjugated_transpose(double complex B[][MAX_SIZE], int rows, int cols, double complex result[][MAX_SIZE]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = B[i][j];
        }
    }
}

// Function for inner product
double inner_product(double a[], double c[], int size) {
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += a[i] * c[i];
    }
    return result;
}

// Function for outer product
void outer_product(double c[], int size_c, double a[], int size_a, double result[][MAX_SIZE]) {
    for (int i = 0; i < size_c; i++) {
        for (int j = 0; j < size_a; j++) {
            result[i][j] = c[i] * a[j];
        }
    }
}

// Function to multiply two matrices
void matrix_mult(double A[][MAX_SIZE], double B[][MAX_SIZE], int rows_A, int cols_A, int cols_B, double result[][MAX_SIZE]) {
    for (int i = 0; i < rows_A; i++) {
        for (int j = 0; j < cols_B; j++) {
            result[i][j] = 0;
            for (int k = 0; k < cols_A; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Function to compute matrix inverse using Gauss-Jordan elimination
void matrix_inverse(double C[][MAX_SIZE], int n, double inv[][MAX_SIZE]) {
    double aug[MAX_SIZE][MAX_SIZE * 2];
    
    // Create augmented matrix [C | I]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = C[i][j];
            aug[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Gauss-Jordan elimination
    for (int i = 0; i < n; i++) {
        // Partial pivoting
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(aug[k][i]) > fabs(aug[max_row][i])) {
                max_row = k;
            }
        }
        
        // Swap rows
        if (max_row != i) {
            for (int j = 0; j < 2 * n; j++) {
                double temp = aug[i][j];
                aug[i][j] = aug[max_row][j];
                aug[max_row][j] = temp;
            }
        }
        
        // Scale pivot row
        double pivot = aug[i][i];
        for (int j = 0; j < 2 * n; j++) {
            aug[i][j] /= pivot;
        }
        
        // Eliminate column
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }
    
    // Extract inverse from augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = aug[i][j + n];
        }
    }
}

// Function to print real matrix
void print_matrix(const char* name, double mat[][MAX_SIZE], int rows, int cols) {
    printf("\n%s:\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.4f ", mat[i][j]);
        }
        printf("\n");
    }
}

// Function to print complex matrix
void print_complex_matrix(const char* name, double complex mat[][MAX_SIZE], int rows, int cols) {
    printf("\n%s:\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("(%6.3f%+6.3fi) ", creal(mat[i][j]), cimag(mat[i][j]));
        }
        printf("\n");
    }
}

int main() {
    srand(time(NULL));
    
    // 1. Matrix Multiplication
    printf("=== Matrix Multiplication ===\n");
    double A[2][MAX_SIZE] = {{1, 2, 3}, {4, 5, 6}};
    double v[3] = {1, 2, 3};
    double product[2];
    
    matrix_vector_mult(A, 2, 3, v, product);
    printf("A * v = [%.2f, %.2f]\n", product[0], product[1]);
    
    // 2. Matrix Creation using Loops
    printf("\n=== Matrix Creation with Random Values ===\n");
    int no_of_rows = 3;
    int no_of_cols = 5;
    double B[MAX_SIZE][MAX_SIZE];
    
    for (int i = 0; i < no_of_rows; i++) {
        for (int j = 0; j < no_of_cols; j++) {
            B[i][j] = (double)rand() / RAND_MAX;
        }
    }
    print_matrix("Matrix B", B, no_of_rows, no_of_cols);
    
    // 3. Sub-matrix Extraction
    printf("\n=== Sub-matrix Extraction ===\n");
    double sub_matrix[2][3];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            sub_matrix[i][j] = B[i][j + 1];  // rows 1-2, cols 2-4 (0-indexed: 0-1, 1-3)
        }
    }
    print_matrix("Sub-matrix of B (rows 1-2, cols 2-4)", sub_matrix, 2, 3);
    
    // 4. Transpose
    printf("\n=== Transpose ===\n");
    double B_transpose[MAX_SIZE][MAX_SIZE];
    transpose(B, no_of_rows, no_of_cols, B_transpose);
    print_matrix("B Transpose", B_transpose, no_of_cols, no_of_rows);
    
    // 5. Complex Matrices
    printf("\n=== Complex Matrices ===\n");
    double complex B_complex[MAX_SIZE][MAX_SIZE];
    for (int i = 0; i < no_of_rows; i++) {
        for (int j = 0; j < no_of_cols; j++) {
            double real_part = (double)rand() / RAND_MAX;
            double imag_part = (double)rand() / RAND_MAX;
            B_complex[i][j] = real_part + I * imag_part;
        }
    }
    print_complex_matrix("Complex Matrix B", B_complex, no_of_rows, no_of_cols);
    
    // 6. Conjugate Transpose
    double complex B_conjugate_transpose[MAX_SIZE][MAX_SIZE];
    conjugate_transpose(B_complex, no_of_rows, no_of_cols, B_conjugate_transpose);
    print_complex_matrix("Conjugate Transpose", B_conjugate_transpose, no_of_cols, no_of_rows);
    
    // 7. Unconjugated Transpose
    double complex B_unconj_transpose[MAX_SIZE][MAX_SIZE];
    unconjugated_transpose(B_complex, no_of_rows, no_of_cols, B_unconj_transpose);
    print_complex_matrix("Unconjugated Transpose", B_unconj_transpose, no_of_cols, no_of_rows);
    
    // 8. Inner and Outer Products
    printf("\n=== Inner and Outer Products ===\n");
    double a[3] = {1, 2, 3};
    double c[3] = {4, 5, 6};
    
    double inner_prod = inner_product(a, c, 3);
    printf("Inner product (a · c) = %.2f\n", inner_prod);
    
    double outer_prod[MAX_SIZE][MAX_SIZE];
    outer_product(c, 3, a, 3, outer_prod);
    print_matrix("Outer product (c ⊗ a)", outer_prod, 3, 3);
    
    // 9. Matrix Inverse
    printf("\n=== Matrix Inverse ===\n");
    double C[MAX_SIZE][MAX_SIZE];
    int n = 8;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = (double)rand() / RAND_MAX;
        }
    }
    
    double inv_C[MAX_SIZE][MAX_SIZE];
    matrix_inverse(C, n, inv_C);
    
    // Verify: C * inv(C) should be identity
    double verify1[MAX_SIZE][MAX_SIZE];
    double verify2[MAX_SIZE][MAX_SIZE];
    matrix_mult(C, inv_C, n, n, n, verify1);
    matrix_mult(inv_C, C, n, n, n, verify2);
    
    print_matrix("C * inv(C) (should be identity)", verify1, n, n);
    print_matrix("inv(C) * C (should be identity)", verify2, n, n);
    
    // Create identity matrix
    double identity[MAX_SIZE][MAX_SIZE];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            identity[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    print_matrix("Identity Matrix (eye(8))", identity, n, n);
    
    return 0;
}
