#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <lapacke.h>

#define MAX_SIZE 100

void print_matrix(const char* name, double* mat, int rows, int cols) {
    printf("\n%s (%dx%d):\n", name, rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.4f ", mat[i * cols + j]);
        }
        printf("\n");
    }
}

int main() {
    srand(time(NULL));
    
    int N = 6;  // rows
    int M = 8;  // columns
    
    printf("=== SVD using LAPACK ===\n");
    
    // Create random matrix A (row-major format)
    double A[N * M];
    for (int i = 0; i < N * M; i++) {
        A[i] = (double)rand() / RAND_MAX;
    }
    
    print_matrix("Original Matrix A", A, N, M);
    
    // Allocate space for SVD results
    double U[N * N];
    double S[N < M ? N : M];  // singular values
    double VT[M * M];
    double superb[N < M ? N - 1 : M - 1];
    
    // Perform SVD
    int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', 
                              N, M, A, M, S, U, N, VT, M, superb);
    
    if (info != 0) {
        printf("SVD failed with error code: %d\n", info);
        return 1;
    }
    
    printf("\n=== SVD Results ===\n");
    print_matrix("U (Left singular vectors)", U, N, N);
    
    printf("\nSingular values:\n");
    for (int i = 0; i < (N < M ? N : M); i++) {
        printf("S[%d] = %.6f\n", i, S[i]);
    }
    
    print_matrix("V' (Right singular vectors transposed)", VT, M, M);
    
    // Reconstruct A = U * Sigma * V'
    double Sigma[N * M];
    for (int i = 0; i < N * M; i++) Sigma[i] = 0.0;
    for (int i = 0; i < (N < M ? N : M); i++) {
        Sigma[i * M + i] = S[i];
    }
    
    double temp[N * M];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            temp[i * M + j] = 0.0;
            for (int k = 0; k < N; k++) {
                temp[i * M + j] += U[i * N + k] * Sigma[k * M + j];
            }
        }
    }
    
    double recreated_A[N * M];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            recreated_A[i * M + j] = 0.0;
            for (int k = 0; k < M; k++) {
                recreated_A[i * M + j] += temp[i * M + k] * VT[k * M + j];
            }
        }
    }
    
    print_matrix("Recreated A from SVD", recreated_A, N, M);
    
    return 0;
}
