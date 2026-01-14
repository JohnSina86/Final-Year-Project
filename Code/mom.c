#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#define MAX_ITER 100000
#define TOL 1e-9

// OPTIMIZED SSOR solver with precomputed inverse diagonals
int ssor_solve_optimized(double complex *Z_flat, double complex *V, double complex *x, 
                         int N, double omega, double tol, int max_iter) {

    double complex *x_new = malloc(N * sizeof(double complex));
    double complex *inv_diag = malloc(N * sizeof(double complex));

    if (!x_new || !inv_diag) {
        fprintf(stderr, "Memory allocation failed\n");
        free(x_new); free(inv_diag);
        return -1;
    }

    // Precompute inverse of diagonal elements
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        inv_diag[i] = 1.0 / Z_flat[i * N + i];
        x[i] = V[i] * inv_diag[i];
        x_new[i] = x[i];
    }

    int converged_iter = -1;

    for (int iter = 0; iter < max_iter; iter++) {
        // Forward sweep - cannot parallelize due to dependencies
        for (int i = 0; i < N; i++) {
            double complex sum = V[i];
            double complex *Z_row = &Z_flat[i * N];

            for (int j = 0; j < i; j++) {
                sum -= Z_row[j] * x_new[j];
            }
            sum -= Z_row[i] * x[i];
            for (int j = i + 1; j < N; j++) {
                sum -= Z_row[j] * x[j];
            }
            x_new[i] = x[i] + omega * sum * inv_diag[i];
        }

        // Backward sweep
        for (int i = N - 1; i >= 0; i--) {
            double complex sum = V[i];
            double complex *Z_row = &Z_flat[i * N];

            for (int j = 0; j < i; j++) {
                sum -= Z_row[j] * x_new[j];
            }
            sum -= Z_row[i] * x_new[i];
            for (int j = i + 1; j < N; j++) {
                sum -= Z_row[j] * x_new[j];
            }
            x_new[i] = x_new[i] + omega * sum * inv_diag[i];
        }

        // Convergence check - parallelized
        double maxdiff = 0.0;
        double maxval = 1e-20;

        #pragma omp parallel for reduction(max:maxdiff,maxval)
        for (int i = 0; i < N; i++) {
            double diff = cabs(x_new[i] - x[i]);
            double val = cabs(x_new[i]);
            if (diff > maxdiff) maxdiff = diff;
            if (val > maxval) maxval = val;
        }

        double relative_error = maxdiff / maxval;

        memcpy(x, x_new, N * sizeof(double complex));

        if (relative_error < tol) {
            converged_iter = iter + 1;
            printf("Converged: %d iterations, error=%.2e\n", converged_iter, relative_error);
            break;
        }

        if ((iter + 1) % 1000 == 0) {
            printf("  Iteration %d: error=%.2e\n", iter + 1, relative_error);
        }
    }

    free(x_new);
    free(inv_diag);
    return converged_iter;
}

// PARALLELIZED OMEGA OPTIMIZATION - Much faster than sequential
double find_omega_parallel(double complex *Z_flat, double complex *V, int N) {
    printf("\n=== Optimizing Relaxation Parameter (Parallel) ===\n");

    // Test 6 omega values in parallel (instead of 20 sequential)
    // NEW (fixed):
    #define NUM_TESTS 6
    double omega_values[NUM_TESTS] = {1.3, 1.4, 1.5, 1.6, 1.7, 1.8};
    double errors[NUM_TESTS];


    printf("Testing %d omega values in parallel...\n", NUM_TESTS);

    // Parallel omega search - each thread tests one omega value
    #pragma omp parallel for schedule(dynamic)
    for (int test = 0; test < NUM_TESTS; test++) {
        double omega = omega_values[test];

        // Allocate per-thread arrays
        double complex *x_test = malloc(N * sizeof(double complex));
        double complex *x_new = malloc(N * sizeof(double complex));
        double complex *inv_diag = malloc(N * sizeof(double complex));

        if (!x_test || !x_new || !inv_diag) {
            free(x_test); free(x_new); free(inv_diag);
            errors[test] = 1e20;
            continue;
        }

        // Precompute inverse diagonals
        for (int i = 0; i < N; i++) {
            inv_diag[i] = 1.0 / Z_flat[i * N + i];
            x_test[i] = V[i] * inv_diag[i];
            x_new[i] = x_test[i];
        }

        // Quick test: 200 iterations (reduced from 500-1000)
        double final_error = 1e20;
        for (int iter = 0; iter < 200; iter++) {
            // Forward sweep
            for (int i = 0; i < N; i++) {
                double complex sum = V[i];
                double complex *Z_row = &Z_flat[i * N];
                for (int j = 0; j < i; j++) sum -= Z_row[j] * x_new[j];
                sum -= Z_row[i] * x_test[i];
                for (int j = i + 1; j < N; j++) sum -= Z_row[j] * x_test[j];
                x_new[i] = x_test[i] + omega * sum * inv_diag[i];
            }

            // Backward sweep
            for (int i = N - 1; i >= 0; i--) {
                double complex sum = V[i];
                double complex *Z_row = &Z_flat[i * N];
                for (int j = 0; j < i; j++) sum -= Z_row[j] * x_new[j];
                sum -= Z_row[i] * x_new[i];
                for (int j = i + 1; j < N; j++) sum -= Z_row[j] * x_new[j];
                x_new[i] = x_new[i] + omega * sum * inv_diag[i];
            }

            // Quick convergence check
            double maxdiff = 0.0, maxval = 1e-20;
            for (int i = 0; i < N; i++) {
                double diff = cabs(x_new[i] - x_test[i]);
                double val = cabs(x_new[i]);
                maxdiff = fmax(maxdiff, diff);
                maxval = fmax(maxval, val);
            }
            final_error = maxdiff / maxval;

            memcpy(x_test, x_new, N * sizeof(double complex));

            // Looser tolerance for quick testing
            if (final_error < 1e-2) break;
        }

        errors[test] = final_error;

        free(x_test);
        free(x_new);
        free(inv_diag);
    }

    // Find best omega from parallel results
    double best_omega = 1.5;
    double best_error = 1e20;

    printf("\nOmega test results:\n");
    for (int i = 0; i < NUM_TESTS; i++) {
        printf("  omega=%.2f: convergence error=%.2e", omega_values[i], errors[i]);
        if (errors[i] < best_error) {
            best_error = errors[i];
            best_omega = omega_values[i];
            printf(" <- BEST");
        }
        printf("\n");
    }

    printf("\nOptimal omega selected: %.2f (test error: %.2e)\n", best_omega, best_error);
    return best_omega;
}

void solve_mom_problem(double length_of_strip, double TX, double TY) {
    clock_t start_time = clock();

    const double MHz = 1000000.0;
    const double MPI = 3.14159265358979323846;
    const double f = 2000.0 * MHz;
    const double epsilon0 = 8.854e-12;
    const double mu0 = 4.0 * MPI * 1.0e-7;
    const double c0 = 1.0 / sqrt(epsilon0 * mu0);
    const double omega = 2.0 * MPI * f;
    const double k0 = omega / c0;
    const double lambda = c0 / f;
    const double eta = sqrt(mu0 / epsilon0);

    /* L-Bracket Geometry */
    double strip1_length = length_of_strip;
    double strip1_angle = 0.0;
    double strip2_length = length_of_strip;
    double strip2_angle = MPI / 2.0;

    int disc_per_lambda = 100;  
    double delta_s = lambda / disc_per_lambda;

    printf("\n==========================================================\n");
    printf("   L-Bracket MoM Analysis (Optimized Implementation)\n");
    printf("==========================================================\n");
    printf("OpenMP threads available: %d\n", omp_get_max_threads());
    printf("Frequency: %.0f MHz\n", f/MHz);
    printf("Wavelength: %.6f m\n", lambda);
    printf("Discretization: %d segments per wavelength\n", disc_per_lambda);
    printf("Segment size: %.6f m (%.4f lambda)\n", delta_s, delta_s/lambda);

    int N_strip1 = (int)ceil(strip1_length / delta_s);
    int N_strip2 = (int)ceil(strip2_length / delta_s);
    int N = N_strip1 + N_strip2;

    printf("\nGeometry:\n");
    printf("  Strip 1 (Horizontal): %d segments\n", N_strip1);
    printf("  Strip 2 (Vertical): %d segments\n", N_strip2);
    printf("  Total segments: N = %d\n", N);
    printf("  Strip length: %.4f m (%.2f lambda)\n", strip1_length, strip1_length/lambda);
    printf("  Transmitter position: (%.4f, %.4f) m\n", TX, TY);

    double matrix_memory_MB = (N * N * sizeof(double complex)) / (1024.0 * 1024.0);
    printf("\nMemory requirements:\n");
    printf("  Impedance matrix: %.2f MB\n", matrix_memory_MB);

    /* Allocate position arrays */
    double *x_pos = malloc(N * sizeof(double));
    double *y_pos = malloc(N * sizeof(double));

    if (!x_pos || !y_pos) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(x_pos); free(y_pos);
        return;
    }

    /* Generate positions - optimized with precomputed values */
    double cos_angle1 = cos(strip1_angle);
    double sin_angle1 = sin(strip1_angle);
    double cos_angle2 = cos(strip2_angle);
    double sin_angle2 = sin(strip2_angle);

    #pragma omp parallel for
    for (int ct = 0; ct < N_strip1; ct++) {
        x_pos[ct] = (ct + 0.5) * cos_angle1 * delta_s;
        y_pos[ct] = (ct + 0.5) * sin_angle1 * delta_s;
    }

    #pragma omp parallel for
    for (int ct = 0; ct < N_strip2; ct++) {
        int idx = N_strip1 + ct;
        x_pos[idx] = (ct + 0.5) * cos_angle2 * delta_s;
        y_pos[idx] = (ct + 0.5) * sin_angle2 * delta_s;
    }

    /* Allocate flat arrays */
    double complex *Z_flat = malloc(N * N * sizeof(double complex));
    double complex *V = malloc(N * sizeof(double complex));
    double complex *x = malloc(N * sizeof(double complex));

    if (!Z_flat || !V || !x) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(x_pos); free(y_pos); free(Z_flat); free(V); free(x);
        return;
    }

    printf("\n=== Building Impedance Matrix (Parallel) ===\n");
    clock_t matrix_start = clock();

    /* Build V vector - parallelized */
    #pragma omp parallel for
    for (int ct = 0; ct < N; ct++) {
        double Rx = x_pos[ct] - TX;
        double Ry = y_pos[ct] - TY;
        double R = sqrt(Rx * Rx + Ry * Ry);
        double J0 = j0(k0 * R);
        double Y0 = y0(k0 * R);
        V[ct] = J0 - I * Y0;
    }

    /* Precompute self term */
    double complex self = -1.0 * delta_s * 
                          (1.0 - I * (2.0 / MPI) * 
                           log(1.781 * k0 * delta_s / (4.0 * exp(1.0))));

    double complex k0_eta_factor = (k0 * eta / 4.0);

    /* Build Z matrix - parallelized */
    printf("Assembling impedance matrix using %d threads...\n", omp_get_max_threads());

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        #pragma omp for schedule(dynamic, 10)
        for (int ct1 = 0; ct1 < N; ct1++) {
            // Progress indicator from thread 0
            if (thread_id == 0 && ct1 % (N / 20) == 0) {
                printf("  Progress: %d%%\r", (int)((ct1 * 100.0) / N));
                fflush(stdout);
            }

            double x1 = x_pos[ct1];
            double y1 = y_pos[ct1];

            for (int ct2 = 0; ct2 < N; ct2++) {
                double complex val;
                if (ct1 != ct2) {
                    double Rx = x1 - x_pos[ct2];
                    double Ry = y1 - y_pos[ct2];
                    double R = sqrt(Rx * Rx + Ry * Ry);
                    double J0 = j0(k0 * R);
                    double Y0 = y0(k0 * R);
                    double complex the_hank = -(J0 - I * Y0);
                    val = k0_eta_factor * the_hank * delta_s;
                } else {
                    val = k0_eta_factor * self;
                }
                Z_flat[ct1 * N + ct2] = val;
            }
        }
    }

    printf("  Progress: 100%% complete\n");
    double matrix_time = (double)(clock() - matrix_start) / CLOCKS_PER_SEC;
    printf("Matrix assembly time: %.2f seconds\n", matrix_time);

    /* Parallel omega optimization */
    clock_t omega_start = clock();
    double optimal_omega = find_omega_parallel(Z_flat, V, N);
    double omega_time = (double)(clock() - omega_start) / CLOCKS_PER_SEC;
    printf("Omega optimization time: %.2f seconds\n", omega_time);

    /* Solve with optimized SSOR */
    printf("\n=== Solving Linear System (SSOR Method) ===\n");
    printf("Using omega = %.4f\n", optimal_omega);
    printf("Convergence tolerance: %.2e\n", TOL);
    clock_t solve_start = clock();

    int iters = ssor_solve_optimized(Z_flat, V, x, N, optimal_omega, TOL, MAX_ITER);

    double solve_time = (double)(clock() - solve_start) / CLOCKS_PER_SEC;

    if (iters == -1) {
        printf("\n*** WARNING: SSOR did not converge within %d iterations ***\n", MAX_ITER);
    } else {
        printf("\n=== Solution Successful ===\n");
        printf("Iterations required: %d\n", iters);
        printf("Solver time: %.2f seconds\n", solve_time);
    }

    double total_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;

    printf("\n==========================================================\n");
    printf("   Computation Summary\n");
    printf("==========================================================\n");
    printf("Matrix assembly:     %.2f seconds\n", matrix_time);
    printf("Omega optimization:  %.2f seconds\n", omega_time);
    printf("SSOR solver:         %.2f seconds\n", solve_time);
    printf("----------------------------------------------------------\n");
    printf("Total time:          %.2f seconds\n", total_time);
    printf("==========================================================\n");

    /* Save results */
    FILE *fx = fopen("current_distribution.txt", "w");
    if (fx) {
        fprintf(fx, "%% ======================================================\n");
        fprintf(fx, "%% Current Distribution Solution\n");
        fprintf(fx, "%% Method: SSOR (Symmetric Successive Over-Relaxation)\n");
        fprintf(fx, "%% ======================================================\n");
        fprintf(fx, "%% Frequency: %.0f MHz\n", f/MHz);
        fprintf(fx, "%% Wavelength: %.6f m\n", lambda);
        fprintf(fx, "%% Resolution: %d segments per wavelength\n", disc_per_lambda);
        fprintf(fx, "%% Geometry: L-Bracket\n");
        fprintf(fx, "%%   N_total=%d, N_horizontal=%d, N_vertical=%d\n", N, N_strip1, N_strip2);
        fprintf(fx, "%% SSOR parameters:\n");
        fprintf(fx, "%%   omega=%.4f (optimized)\n", optimal_omega);
        fprintf(fx, "%%   iterations=%d\n", iters);
        fprintf(fx, "%%   tolerance=%.2e\n", TOL);
        fprintf(fx, "%% Computation time: %.2f seconds\n", total_time);
        fprintf(fx, "%% ======================================================\n");
        fprintf(fx, "%% Format: (real, imaginary) - current coefficients\n");
        fprintf(fx, "%% ======================================================\n");
        for (int i = 0; i < N; i++) {
            fprintf(fx, "(%15.8e,%15.8e)\n", creal(x[i]), cimag(x[i]));
        }
        fclose(fx);
        printf("\nResults saved to: current_distribution.txt\n");
    }

    FILE *fg = fopen("geometry.txt", "w");
    if (fg) {
        fprintf(fg, "%% Segment center positions (x, y) in meters\n");
        fprintf(fg, "%% Resolution: %d segments per wavelength\n", disc_per_lambda);
        fprintf(fg, "%% N_total=%d, N_horizontal=%d, N_vertical=%d\n", N, N_strip1, N_strip2);
        fprintf(fg, "%% Segments 0-%d: Horizontal strip (along x-axis)\n", N_strip1-1);
        fprintf(fg, "%% Segments %d-%d: Vertical strip (along y-axis)\n", N_strip1, N-1);
        for (int i = 0; i < N; i++) {
            fprintf(fg, "%.8e %.8e\n", x_pos[i], y_pos[i]);
        }
        fclose(fg);
        printf("Geometry saved to: geometry.txt\n");
    }

    /* Cleanup */
    free(x_pos); free(y_pos);
    free(Z_flat); free(V); free(x);
}

int main(void) {
    double length_of_strip, TX, TY;

    FILE *f_pos = fopen("Transmitter pos", "r");
    if (f_pos == NULL) {
        fprintf(stderr, "Error: Cannot open 'Transmitter pos' file.\n");
        fprintf(stderr, "Please create this file with:\n");
        fprintf(stderr, "  Line 1: Strip length (m)\n");
        fprintf(stderr, "  Line 2: Transmitter X position (m)\n");
        fprintf(stderr, "  Line 3: Transmitter Y position (m)\n");
        return 1;
    }

    if (fscanf(f_pos, "%lf", &length_of_strip) != 1 ||
        fscanf(f_pos, "%lf", &TX) != 1 ||
        fscanf(f_pos, "%lf", &TY) != 1) {
        fprintf(stderr, "Error: Invalid format in 'Transmitter pos' file.\n");
        fclose(f_pos);
        return 1;
    }
    fclose(f_pos);

    solve_mom_problem(length_of_strip, TX, TY);

    printf("\n*** Computation complete ***\n\n");
    return 0;
}
