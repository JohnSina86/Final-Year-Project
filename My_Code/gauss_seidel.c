#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 100

// Function to perform Gauss-Seidel iteration
void gauss_seidel(double A[][MAX_SIZE], double b[], double x[], int N, int iterations) {
    double sum;
    
    // Initialize x to zeros
    for (int i = 0; i < N; i++) {
        x[i] = 0.0;
    }
    
    // Iterate for the specified number of times
    for (int iter = 0; iter < iterations; iter++) {
        for (int i = 0; i < N; i++) {
            sum = b[i];
            
            // Subtract contributions from already updated components (current iteration)
            for (int j = 0; j < i; j++) {
                sum -= A[i][j] * x[j];
            }
            
            // Subtract contributions from previous iteration components
            for (int j = i + 1; j < N; j++) {
                sum -= A[i][j] * x[j];
            }
            
            // Update x[i]
            x[i] = sum / A[i][i];
        }
    }
}

// Function to calculate exact solution using Gaussian elimination
void gaussian_elimination(double A[][MAX_SIZE], double b[], double x[], int N) {
    double temp[MAX_SIZE][MAX_SIZE];
    double b_temp[MAX_SIZE];
    
    // Copy A and b to temporary arrays
    for (int i = 0; i < N; i++) {
        b_temp[i] = b[i];
        for (int j = 0; j < N; j++) {
            temp[i][j] = A[i][j];
        }
    }
    
    // Forward elimination
    for (int i = 0; i < N; i++) {
        // Partial pivoting
        int max_row = i;
        for (int k = i + 1; k < N; k++) {
            if (fabs(temp[k][i]) > fabs(temp[max_row][i])) {
                max_row = k;
            }
        }
        
        // Swap rows
        if (max_row != i) {
            for (int j = 0; j < N; j++) {
                double t = temp[i][j];
                temp[i][j] = temp[max_row][j];
                temp[max_row][j] = t;
            }
            double t = b_temp[i];
            b_temp[i] = b_temp[max_row];
            b_temp[max_row] = t;
        }
        
        // Eliminate column
        for (int k = i + 1; k < N; k++) {
            double factor = temp[k][i] / temp[i][i];
            for (int j = i; j < N; j++) {
                temp[k][j] -= factor * temp[i][j];
            }
            b_temp[k] -= factor * b_temp[i];
        }
    }
    
    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
        x[i] = b_temp[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= temp[i][j] * x[j];
        }
        x[i] /= temp[i][i];
    }
}

// Function to print vector
void print_vector(const char* name, double v[], int N) {
    printf("\n%s:\n", name);
    for (int i = 0; i < N; i++) {
        printf("%10.6f ", fabs(v[i]));
        if ((i + 1) % 5 == 0) printf("\n");
    }
    printf("\n");
}

int main() {
    int N = 4;  // Size of the system
    int no_of_iterations = 10;
    
    // Example system: 4x + y = 5, x + 3y = 6, 2x + y + 5z = 8, x + z + 4w = 7
    double A[MAX_SIZE][MAX_SIZE] = {
        {4, 1, 0, 0},
        {1, 3, 0, 0},
        {2, 1, 5, 0},
        {1, 0, 1, 4}
    };
    
    double b[MAX_SIZE] = {5, 6, 8, 7};
    
    double x_GS[MAX_SIZE];
    double exact_x[MAX_SIZE];
    
    // Calculate exact solution
    gaussian_elimination(A, b, exact_x, N);
    
    // Calculate Gauss-Seidel solution
    gauss_seidel(A, b, x_GS, N, no_of_iterations);
    
    // Print results
    printf("Number of Gauss-Seidel iterations: %d\n", no_of_iterations);
    print_vector("Gauss-Seidel Solution (absolute values)", x_GS, N);
    print_vector("Exact Solution (absolute values)", exact_x, N);
    
    // Calculate and print error
    double error = 0.0;
    for (int i = 0; i < N; i++) {
        error += fabs(x_GS[i] - exact_x[i]);
    }
    printf("\nTotal absolute error: %.6e\n", error);
    
    return 0;
}
