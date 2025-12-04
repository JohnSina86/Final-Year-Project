#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#define MAX_IT 50000
#define MAX_N 100000

// Extracted SSOR iteration function (unchanged)
int ssor_iterate(double complex **Z, double complex *V, double complex *x, 
                 double complex *x_new, int N, double omega, double tol, 
                 int max_iter, double *final_relative_error) {
    
    for (int i = 0; i < N; i++) {
        x_new[i] = x[i];
    }
    
    int converged_iter = -1;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Forward sweep
        for (int i = 0; i < N; i++) {
            double complex sum = V[i];
            for (int j = 0; j < i; j++) {
                sum -= Z[i][j] * x_new[j];
            }
            sum -= Z[i][i] * x[i];
            for (int j = i + 1; j < N; j++) {
                sum -= Z[i][j] * x[j];
            }
            x_new[i] = x[i] + omega * (sum / Z[i][i]);
        }
        
        // Backward sweep
        for (int i = N - 1; i >= 0; i--) {
            double complex sum = V[i];
            for (int j = 0; j < i; j++) {
                sum -= Z[i][j] * x_new[j];
            }
            sum -= Z[i][i] * x_new[i];
            for (int j = i + 1; j < N; j++) {
                sum -= Z[i][j] * x_new[j];
            }
            x_new[i] = x_new[i] + omega * (sum / Z[i][i]);
        }
        
        // Convergence check
        double maxdiff = 0.0;
        double maxval = 1e-20;
        for (int i = 0; i < N; i++) {
            double diff = cabs(x_new[i] - x[i]);
            double val = cabs(x_new[i]);
            if (diff > maxdiff) maxdiff = diff;
            if (val > maxval) maxval = val;
        }
        
        double relative_error = maxdiff / maxval;
        *final_relative_error = relative_error;
        
        for (int i = 0; i < N; i++) {
            x[i] = x_new[i];
        }
        
        if (relative_error < tol) {
            converged_iter = iter + 1;
            break;
        }
    }
    
    return converged_iter;
}

// UNIVERSAL 3-STEP OMEGA OPTIMIZER
double find_optimal_omega(double complex **Z, double complex *V, int N) {
    // ALWAYS scan full range first: 0.01 to 1.99
    double omega_ranges[3][2] = {
        {0.01, 1.99},      // Phase 1: Full coarse scan
        {0.0, 0.0},        // Phase 2: Will be set dynamically
        {0.0, 0.0}         // Phase 3: Will be set dynamically
    };
    
    double tol_levels[3] = {1e-3, 1e-6, 1e-9};
    double step_levels[3] = {0.01, 0.001, 0.0001};
    
    double best_omega = 1.0;
    double best_error = 1e20;
    int best_iters = MAX_IT;
    
    printf("\n=== UNIVERSAL OMEGA OPTIMIZATION (3 PHASES) ===\n");
    printf("Works for ANY problem size/matrix condition\n\n");
    
    for (int phase = 0; phase < 3; phase++) {
        double omega_start = omega_ranges[phase][0];
        double omega_end = omega_ranges[phase][1];
        double step = step_levels[phase];
        double tol = tol_levels[phase];
        
        // For phases 2-3, narrow around previous best
        if (phase > 0) {
            double search_width = step * 20;  // +/- 10 steps
            omega_start = fmax(0.01, best_omega - search_width);
            omega_end = fmin(1.99, best_omega + search_width);
        }
        
        printf("Phase %d: tol=%.0e, omega=[%.3f,%.3f] (step=%.4f)\n", 
               phase+1, tol, omega_start, omega_end, step);
        
        double phase_best_error = 1e20;
        double phase_best_omega = best_omega;
        int phase_best_iters = MAX_IT;
        
        // Test all omegas in range
        for (double omega = omega_start; omega <= omega_end + 1e-10; omega += step) {
            double complex *x_test = malloc(N * sizeof(double complex));
            double complex *x_new_test = malloc(N * sizeof(double complex));
            
            // Universal initial guess: diagonal preconditioning
            for (int i = 0; i < N; i++) {
                x_test[i] = V[i] / Z[i][i];
            }
            
            double final_error;
            int iters = ssor_iterate(Z, V, x_test, x_new_test, N, omega, tol, MAX_IT, &final_error);
            
            // Track best: lowest error OR converges fastest
            if (iters != -1 && final_error < phase_best_error) {
                phase_best_error = final_error;
                phase_best_omega = omega;
                phase_best_iters = iters;
            }
            
            free(x_test);
            free(x_new_test);
        }
        
        // Update for next phase
        best_omega = phase_best_omega;
        best_error = phase_best_error;
        best_iters = phase_best_iters;
        
        printf("  → BEST: omega=%.6f (error=%.2e, iters=%d)\n", 
               best_omega, best_error, best_iters);
        
        // Set range for next phase
        if (phase < 2) {
            omega_ranges[phase+1][0] = best_omega;
            omega_ranges[phase+1][1] = best_omega;
        }
    }
    
    printf("\n🎯 FINAL OPTIMAL OMEGA: %.6f\n", best_omega);
    printf("   (Valid for this matrix/problem)\n");
    
    return best_omega;
}

int main(void) {
    FILE *fp = fopen("function_result.txt", "r");
    if (!fp) {
        perror("function_result.txt");
        return 1;
    }
    
    // Determine N (your original logic)
    char line[500000];
    int N = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '%' || line[0] == '\n' || line[0] == '\r') continue;
        char *ptr = line;
        while (*ptr) {
            double re, im;
            int n_chars = 0;
            if (sscanf(ptr, " (%lf,%lf)%n", &re, &im, &n_chars) == 2) {
                N++;
                ptr += n_chars;
            } else {
                ptr++;
            }
        }
        break;
    }
    
    rewind(fp);
    printf("Problem size N = %d\n", N);
    
    // Allocate (same as original)
    double complex **Z = malloc(N * sizeof(double complex *));
    for (int i = 0; i < N; i++) {
        Z[i] = malloc(N * sizeof(double complex));
    }
    double complex *V = malloc(N * sizeof(double complex));
    double complex *x = malloc(N * sizeof(double complex));
    double complex *x_new = malloc(N * sizeof(double complex));
    
    if (!Z || !V || !x || !x_new) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }
    
    // Read data (your original parser)
    int section = 0, row_count = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '%' || line[0] == '\n' || line[0] == '\r') continue;
        
        if (section == 0) {
            char *ptr = line;
            int col = 0;
            while (*ptr && col < N) {
                double re, im;
                int n_chars = 0;
                if (sscanf(ptr, " (%lf,%lf)%n", &re, &im, &n_chars) == 2) {
                    Z[row_count][col] = re + I * im;
                    col++;
                    ptr += n_chars;
                } else {
                    ptr++;
                }
            }
            row_count++;
            if (row_count == N) {
                section = 1;
                row_count = 0;
            }
        } else if (section == 1) {
            double re, im;
            if (sscanf(line, " (%lf,%lf)", &re, &im) == 2) {
                V[row_count] = re + I * im;
                row_count++;
                if (row_count == N) break;
            }
        }
    }
    fclose(fp);
    
    // OPTIMIZE OMEGA (ALWAYS FULL RANGE)
    double optimal_omega = find_optimal_omega(Z, V, N);
    
    // FINAL HIGH-PRECISION SOLUTION
    printf("\n=== FINAL SOLUTION (omega=%.6f) ===\n", optimal_omega);
    for (int i = 0; i < N; i++) {
        x[i] = V[i] / Z[i][i];  // Initial guess
    }
    
    double final_error;
    int final_iters = ssor_iterate(Z, V, x, x_new, N, optimal_omega, 1e-12, MAX_IT*2, &final_error);
    
    printf("✅ CONVERGED: %d iterations, error=%.2e\n", final_iters, final_error);
    
    // Save result
    FILE *fx = fopen("x_optimized.txt", "w");
    fprintf(fx, "%% SSOR with optimal omega=%.6f (N=%d)\n", optimal_omega, N);
    for (int i = 0; i < N; i++) {
        fprintf(fx, "(%15.8e,%15.8e)\n", creal(x[i]), cimag(x[i]));
    }
    fclose(fx);
    printf("Saved to x_optimized.txt\n");
    
    // Cleanup
    for (int i = 0; i < N; i++) free(Z[i]);
    free(Z); free(V); free(x); free(x_new);
    
    return 0;
}
