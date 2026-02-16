#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <malloc.h>

/*******************************************************************************
 * L-BRACKET METHOD OF MOMENTS (MoM) SOLVER
 * 
 * Optimized SSOR iterative solver with hierarchical omega optimization.
 * Features:
 *   - Parallel matrix assembly with cache blocking
 *   - 3-round hierarchical omega search (coarse → fine → ultra-fine)
 *   - Parallel SSOR solver with convergence detection
 *   - Optimized for 8-core systems
 * 
 * Compile: gcc -std=c11 -O3 -march=native -fopenmp -ffast-math -funroll-loops \
 *              mom_clean.c -lm -o mom_clean.exe
 ******************************************************************************/

//==============================================================================
// CONFIGURATION CONSTANTS
//==============================================================================

// Solver parameters
#define MAX_ITERATIONS      100000      // Maximum SSOR iterations
#define CONVERGENCE_TOL     1e-9        // Relative error threshold
#define CACHE_BLOCK_SIZE    64          // Cache-friendly block size

// Omega search parameters
#define OMEGA_MIN           0.1         // Minimum relaxation parameter
#define OMEGA_MAX           1.9         // Maximum relaxation parameter
#define OMEGA_TEST_ITERS    1000         // Iterations per omega test

// Physical constants
#define M_PI 3.14159265358979323846
#define FREQUENCY_MHZ       2400.0      // Operating frequency (MHz)
#define EPSILON_0           8.854e-12   // Vacuum permittivity (F/m)
#define MU_0                1.256637e-6 // Vacuum permeability (H/m)
#define SEGMENTS_PER_LAMBDA 100         // Discretization density

//==============================================================================
// PLATFORM COMPATIBILITY
//==============================================================================

#ifdef _WIN32
    #define aligned_alloc(alignment, size) _aligned_malloc(size, alignment)
    #define aligned_free(ptr) _aligned_free(ptr)
    #define bessel_j0(x) _j0(x)
    #define bessel_y0(x) _y0(x)
#else
    #define aligned_free(ptr) free(ptr)
    #define bessel_j0(x) j0(x)
    #define bessel_y0(x) y0(x)
#endif

//==============================================================================
// DATA STRUCTURES
//==============================================================================

/**
 * Electromagnetic problem parameters
 */
typedef struct {
    double frequency;           // Hz
    double wavelength;          // m
    double wavenumber;          // rad/m
    double impedance;           // ohms (free space)
    double segment_length;      // m
} EMParameters;

/**
 * Geometry configuration
 */
typedef struct {
    double strip_length;        // Length of each strip (m)
    int num_horizontal;         // Segments in horizontal strip
    int num_vertical;           // Segments in vertical strip
    int total_segments;         // Total N
    double *x_positions;        // Segment centers X (m)
    double *y_positions;        // Segment centers Y (m)
} Geometry;

/**
 * Transmitter location
 */
typedef struct {
    double x;                   // m
    double y;                   // m
} Transmitter;

/**
 * Solver convergence result
 */
typedef struct {
    int iterations;
    double final_error;
    double time_seconds;
    int converged;
} ConvergenceResult;

//==============================================================================
// ELECTROMAGNETIC CALCULATIONS
//==============================================================================

/**
 * Initialize electromagnetic parameters based on frequency
 */
EMParameters initialize_em_parameters(double frequency_mhz) {
    EMParameters params;

    params.frequency = frequency_mhz * 1e6;  // Convert to Hz

    double speed_of_light = 1.0 / sqrt(EPSILON_0 * MU_0);
    params.wavelength = speed_of_light / params.frequency;
    params.wavenumber = 2.0 * M_PI / params.wavelength;
    params.impedance = sqrt(MU_0 / EPSILON_0);
    params.segment_length = params.wavelength / SEGMENTS_PER_LAMBDA;

    return params;
}

/**
 * Calculate self-impedance element for a segment
 */
double complex calculate_self_impedance(EMParameters params) {
    double delta_s = params.segment_length;
    double k0 = params.wavenumber;

    // Analytical self-impedance formula for thin wire
    double complex self_term = -delta_s * (1.0 - I * (2.0 / M_PI) * 
                               log(1.781 * k0 * delta_s / (4.0 * exp(1.0))));

    return (params.wavenumber * params.impedance / 4.0) * self_term;
}

/**
 * Calculate mutual impedance between two segments
 */
double complex calculate_mutual_impedance(double x1, double y1, double x2, double y2,
                                         EMParameters params) {
    double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    double k_times_r = params.wavenumber * distance;

    // Hankel function approximation: H0(2) ≈ -(J0 - i*Y0)
    double J0 = bessel_j0(k_times_r);
    double Y0 = bessel_y0(k_times_r);
    double complex hankel = -(J0 - I * Y0);

    return (params.wavenumber * params.impedance / 4.0) * 
           hankel * params.segment_length;
}

/**
 * Calculate excitation voltage at a segment due to transmitter
 */
double complex calculate_excitation(double seg_x, double seg_y, 
                                   Transmitter tx, EMParameters params) {
    double distance = sqrt((seg_x - tx.x) * (seg_x - tx.x) + 
                          (seg_y - tx.y) * (seg_y - tx.y));
    double k_times_r = params.wavenumber * distance;

    double J0 = bessel_j0(k_times_r);
    double Y0 = bessel_y0(k_times_r);

    return J0 - I * Y0;
}

//==============================================================================
// GEOMETRY SETUP
//==============================================================================

/**
 * Create L-bracket geometry with horizontal and vertical strips
 */
Geometry create_geometry(double strip_length, EMParameters params) {
    Geometry geom;

    geom.strip_length = strip_length;
    geom.num_horizontal = (int)ceil(strip_length / params.segment_length);
    geom.num_vertical = (int)ceil(strip_length / params.segment_length);
    geom.total_segments = geom.num_horizontal + geom.num_vertical;

    // Allocate position arrays
    geom.x_positions = aligned_alloc(64, geom.total_segments * sizeof(double));
    geom.y_positions = aligned_alloc(64, geom.total_segments * sizeof(double));

    // Horizontal strip (along x-axis)
    #pragma omp parallel for
    for (int i = 0; i < geom.num_horizontal; i++) {
        geom.x_positions[i] = (i + 0.5) * params.segment_length;
        geom.y_positions[i] = 0.0;
    }

    // Vertical strip (along y-axis)
    #pragma omp parallel for
    for (int i = 0; i < geom.num_vertical; i++) {
        int idx = geom.num_horizontal + i;
        geom.x_positions[idx] = 0.0;
        geom.y_positions[idx] = (i + 0.5) * params.segment_length;
    }

    return geom;
}

/**
 * Free geometry memory
 */
void free_geometry(Geometry *geom) {
    aligned_free(geom->x_positions);
    aligned_free(geom->y_positions);
}

//==============================================================================
// MATRIX ASSEMBLY
//==============================================================================

/**
 * Build impedance matrix Z and excitation vector V
 * Uses cache blocking and parallel assembly for performance
 */
void assemble_system(double complex *Z, double complex *V, 
                    Geometry geom, Transmitter tx, EMParameters params) {

    int N = geom.total_segments;
    double complex self_impedance = calculate_self_impedance(params);

    printf("\nAssembling %dx%d impedance matrix using %d threads...\n", 
           N, N, omp_get_max_threads());

    // Build excitation vector (parallel transmitter interaction)
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        V[i] = calculate_excitation(geom.x_positions[i], geom.y_positions[i], 
                                    tx, params);
    }

    // Build impedance matrix with cache blocking
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        #pragma omp for schedule(dynamic)
        for (int i_block = 0; i_block < N; i_block += CACHE_BLOCK_SIZE) {
            int i_end = fmin(i_block + CACHE_BLOCK_SIZE, N);

            // Progress indicator (thread 0 only)
            if (thread_id == 0 && i_block % (N / 20 + 1) == 0) {
                printf("  Progress: %d%%\r", (int)(i_block * 100.0 / N));
                fflush(stdout);
            }

            // Process block of rows
            for (int i = i_block; i < i_end; i++) {
                double xi = geom.x_positions[i];
                double yi = geom.y_positions[i];

                // Process in column blocks for cache efficiency
                for (int j_block = 0; j_block < N; j_block += CACHE_BLOCK_SIZE) {
                    int j_end = fmin(j_block + CACHE_BLOCK_SIZE, N);

                    for (int j = j_block; j < j_end; j++) {
                        if (i == j) {
                            Z[i * N + j] = self_impedance;
                        } else {
                            Z[i * N + j] = calculate_mutual_impedance(
                                xi, yi, geom.x_positions[j], geom.y_positions[j], params);
                        }
                    }
                }
            }
        }
    }

    printf("  Progress: 100%% complete\n");
}

//==============================================================================
// SSOR ITERATION CORE
//==============================================================================

/**
 * Perform one SSOR iteration (forward + backward sweep)
 * Returns: relative error after this iteration
 */
double ssor_iteration(double complex *Z, double complex *V, 
                     double complex *x_old, double complex *x_new,
                     double complex *inv_diag, double complex *diag_vals,
                     int N, double omega) {

    // FORWARD SWEEP (i = 0 to N-1)
    for (int i = 0; i < N; i++) {
        double complex sum = V[i];
        double complex *Z_row = &Z[i * N];

        // Lower triangle (uses new values)
        for (int j = 0; j < i; j++) {
            sum -= Z_row[j] * x_new[j];
        }

        // Diagonal (uses old value)
        sum -= diag_vals[i] * x_old[i];

        // Upper triangle (uses old values)
        for (int j = i + 1; j < N; j++) {
            sum -= Z_row[j] * x_old[j];
        }

        // Update with relaxation
        x_new[i] = x_old[i] + omega * sum * inv_diag[i];
    }

    // BACKWARD SWEEP (i = N-1 to 0)
    for (int i = N - 1; i >= 0; i--) {
        double complex sum = V[i];
        double complex *Z_row = &Z[i * N];

        // Lower triangle
        for (int j = 0; j < i; j++) {
            sum -= Z_row[j] * x_new[j];
        }

        // Diagonal (uses current value)
        sum -= diag_vals[i] * x_new[i];

        // Upper triangle
        for (int j = i + 1; j < N; j++) {
            sum -= Z_row[j] * x_new[j];
        }

        // Second relaxation step
        x_new[i] = x_new[i] + omega * sum * inv_diag[i];
    }

    // Calculate convergence metric (parallel reduction)
    double max_change = 0.0;
    double max_magnitude = 1e-20;

    #pragma omp parallel
    {
        double local_max_change = 0.0;
        double local_max_mag = 1e-20;

        #pragma omp for nowait
        for (int i = 0; i < N; i++) {
            double change = cabs(x_new[i] - x_old[i]);
            double magnitude = cabs(x_new[i]);
            local_max_change = fmax(local_max_change, change);
            local_max_mag = fmax(local_max_mag, magnitude);
        }

        #pragma omp critical
        {
            max_change = fmax(max_change, local_max_change);
            max_magnitude = fmax(max_magnitude, local_max_mag);
        }
    }

    return max_change / max_magnitude;
}

//==============================================================================
// OMEGA OPTIMIZATION
//==============================================================================

/**
 * Test convergence rate for a single omega value
 * Returns final error after test_iterations
 */
double test_omega_value(double complex *Z, double complex *V, int N,
                       double omega, int test_iterations) {

    // Allocate temporary arrays
    double complex *x = aligned_alloc(64, N * sizeof(double complex));
    double complex *x_new = aligned_alloc(64, N * sizeof(double complex));
    double complex *inv_diag = aligned_alloc(64, N * sizeof(double complex));
    double complex *diag = aligned_alloc(64, N * sizeof(double complex));

    if (!x || !x_new || !inv_diag || !diag) {
        aligned_free(x); aligned_free(x_new);
        aligned_free(inv_diag); aligned_free(diag);
        return 1e20;  // Failed allocation
    }

    // Initialize: x0 = D^(-1) * V (diagonal preconditioner)
    for (int i = 0; i < N; i++) {
        diag[i] = Z[i * N + i];
        inv_diag[i] = 1.0 / diag[i];
        x[i] = V[i] * inv_diag[i];
        x_new[i] = x[i];
    }

    // Run fixed number of iterations
    double final_error = 1e20;
    for (int iter = 0; iter < test_iterations; iter++) {
        double error = ssor_iteration(Z, V, x, x_new, inv_diag, diag, N, omega);

        // Sample error periodically
        if ((iter + 1) % 100 == 0 || iter == test_iterations - 1) {
            final_error = error;
        }

        // Copy x_new → x for next iteration
        memcpy(x, x_new, N * sizeof(double complex));
    }

    aligned_free(x); aligned_free(x_new);
    aligned_free(inv_diag); aligned_free(diag);

    return final_error;
}

/**
 * Find best omega in a range using parallel testing
 * Returns index of best omega
 */
int find_best_omega(double *omega_values, double *errors, int count,
                   double complex *Z, double complex *V, int N, int test_iters) {

    // Test all omega values in parallel
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < count; i++) {
        errors[i] = test_omega_value(Z, V, N, omega_values[i], test_iters);
    }

    // Find minimum error
    int best_idx = 0;
    for (int i = 1; i < count; i++) {
        if (errors[i] < errors[best_idx]) {
            best_idx = i;
        }
    }

    return best_idx;
}

/**
 * Three-round hierarchical omega search
 * Round 1: Coarse (0.1 steps)
 * Round 2: Fine (0.01 steps)
 * Round 3: Ultra-fine (0.001 steps)
 */
double optimize_omega_hierarchical(double complex *Z, double complex *V, int N) {
    printf("\n=== Hierarchical Omega Optimization ===\n");

    // ROUND 1: Coarse search [1.0, 1.1, ..., 1.9]
    printf("\nRound 1: Coarse search (0.1 increments)\n");

    const int round1_count = 19;
    double omega_r1[19], errors_r1[19];

    for (int i = 0; i < round1_count; i++) {
        omega_r1[i] = OMEGA_MIN + i * 0.1;
    }

    int best_r1 = find_best_omega(omega_r1, errors_r1, round1_count, 
                                   Z, V, N, OMEGA_TEST_ITERS);

    printf("  Best: omega=%.2f (error=%.4e)\n", omega_r1[best_r1], errors_r1[best_r1]);

    // Select neighbor for Round 2 range
    int neighbor_r1 = (best_r1 == 0) ? 1 : 
                      (best_r1 == round1_count - 1) ? round1_count - 2 :
                      (errors_r1[best_r1 - 1] < errors_r1[best_r1 + 1]) ? best_r1 - 1 : best_r1 + 1;

    double omega_min_r2 = fmin(omega_r1[best_r1], omega_r1[neighbor_r1]);
    double omega_max_r2 = fmax(omega_r1[best_r1], omega_r1[neighbor_r1]);

    // ROUND 2: Fine search
    printf("\nRound 2: Fine search [%.2f to %.2f] (0.01 increments)\n", 
           omega_min_r2, omega_max_r2);

    const int round2_count = 11;
    double omega_r2[11], errors_r2[11];

    for (int i = 0; i < round2_count; i++) {
        omega_r2[i] = omega_min_r2 + i * 0.01;
    }

    int best_r2 = find_best_omega(omega_r2, errors_r2, round2_count,
                                   Z, V, N, OMEGA_TEST_ITERS);

    printf("  Best: omega=%.3f (error=%.4e)\n", omega_r2[best_r2], errors_r2[best_r2]);

    // Select neighbor for Round 3
    int neighbor_r2 = (best_r2 == 0) ? 1 :
                      (best_r2 == round2_count - 1) ? round2_count - 2 :
                      (errors_r2[best_r2 - 1] < errors_r2[best_r2 + 1]) ? best_r2 - 1 : best_r2 + 1;

    double omega_min_r3 = fmin(omega_r2[best_r2], omega_r2[neighbor_r2]);
    double omega_max_r3 = fmax(omega_r2[best_r2], omega_r2[neighbor_r2]);

    // ROUND 3: Ultra-fine search
    printf("\nRound 3: Ultra-fine search [%.3f to %.3f] (0.001 increments)\n",
           omega_min_r3, omega_max_r3);

    const int round3_count = 11;
    double omega_r3[11], errors_r3[11];

    for (int i = 0; i < round3_count; i++) {
        omega_r3[i] = omega_min_r3 + i * 0.001;
    }

    int best_r3 = find_best_omega(omega_r3, errors_r3, round3_count,
                                   Z, V, N, OMEGA_TEST_ITERS);

    double optimal_omega = omega_r3[best_r3];

    printf("\n==========================================================\n");
    printf("  OPTIMAL OMEGA: %.4f (error: %.4e)\n", optimal_omega, errors_r3[best_r3]);
    printf("==========================================================\n");

    return optimal_omega;
}

//==============================================================================
// MAIN SOLVER
//==============================================================================

/**
 * Solve the system using SSOR with given omega
 */
ConvergenceResult solve_ssor(double complex *Z, double complex *V, double complex *x,
                             int N, double omega, double tolerance, int max_iter) {

    ConvergenceResult result = {0, 1e20, 0, 0};
    clock_t start = clock();

    // Allocate working arrays
    double complex *x_new = aligned_alloc(64, N * sizeof(double complex));
    double complex *inv_diag = aligned_alloc(64, N * sizeof(double complex));
    double complex *diag = aligned_alloc(64, N * sizeof(double complex));

    if (!x_new || !inv_diag || !diag) {
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        aligned_free(x_new); aligned_free(inv_diag); aligned_free(diag);
        return result;
    }

    // Initialize
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        diag[i] = Z[i * N + i];
        inv_diag[i] = 1.0 / diag[i];
        x[i] = V[i] * inv_diag[i];
        x_new[i] = x[i];
    }

    // Iteration loop
    printf("\nSolving with SSOR (omega=%.4f, tol=%.2e)...\n", omega, tolerance);

    for (int iter = 0; iter < max_iter; iter++) {
        double error = ssor_iteration(Z, V, x, x_new, inv_diag, diag, N, omega);
        result.final_error = error;

        // Check convergence
        if (error < tolerance) {
            result.iterations = iter + 1;
            result.converged = 1;
            printf("Converged in %d iterations (error=%.2e)\n", 
                   result.iterations, error);
            break;
        }

        // Progress report
        if ((iter + 1) % 1000 == 0) {
            printf("  Iteration %d: error=%.2e\n", iter + 1, error);
        }

        // Update for next iteration
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            x[i] = x_new[i];
        }
    }

    if (!result.converged) {
        printf("WARNING: Did not converge in %d iterations\n", max_iter);
    }

    result.time_seconds = (double)(clock() - start) / CLOCKS_PER_SEC;

    aligned_free(x_new); aligned_free(inv_diag); aligned_free(diag);
    return result;
}

//==============================================================================
// FILE I/O
//==============================================================================

/**
 * Save solution to file
 */
void save_results(const char *filename, double complex *x, int N,
                 EMParameters params, Geometry geom, Transmitter tx,
                 double omega, ConvergenceResult result, double total_time) {

    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "WARNING: Could not save results to %s\n", filename);
        return;
    }

    fprintf(f, "%% ======================================================\n");
    fprintf(f, "%% L-Bracket Method of Moments Solution\n");
    fprintf(f, "%% ======================================================\n");
    fprintf(f, "%% Frequency: %.0f MHz\n", params.frequency / 1e6);
    fprintf(f, "%% Wavelength: %.6f m\n", params.wavelength);
    fprintf(f, "%% Segments: N=%d (%d horiz + %d vert)\n", 
            N, geom.num_horizontal, geom.num_vertical);
    fprintf(f, "%% Transmitter: (%.4f, %.4f) m\n", tx.x, tx.y);
    fprintf(f, "%%\n");
    fprintf(f, "%% Solver: SSOR with hierarchical omega\n");
    fprintf(f, "%%   omega = %.4f\n", omega);
    fprintf(f, "%%   iterations = %d\n", result.iterations);
    fprintf(f, "%%   convergence = %.2e\n", result.final_error);
    fprintf(f, "%%   time = %.2f seconds\n", total_time);
    fprintf(f, "%% ======================================================\n");
    fprintf(f, "%% Format: (real, imaginary) current coefficients\n");
    fprintf(f, "%% ======================================================\n");

    for (int i = 0; i < N; i++) {
        fprintf(f, "(%15.8e,%15.8e)\n", creal(x[i]), cimag(x[i]));
    }

    fclose(f);
}

/**
 * Save geometry to file
 */
void save_geometry(const char *filename, Geometry geom, EMParameters params) {
    FILE *f = fopen(filename, "w");
    if (!f) return;

    fprintf(f, "%% L-Bracket Geometry\n");
    fprintf(f, "%% N_total=%d, N_horizontal=%d, N_vertical=%d\n",
            geom.total_segments, geom.num_horizontal, geom.num_vertical);
    fprintf(f, "%% Resolution: %d segments per wavelength\n", SEGMENTS_PER_LAMBDA);
    fprintf(f, "%% Format: x(m) y(m)\n");

    for (int i = 0; i < geom.total_segments; i++) {
        fprintf(f, "%.8e %.8e\n", geom.x_positions[i], geom.y_positions[i]);
    }

    fclose(f);
}

//==============================================================================
// MAIN PROGRAM
//==============================================================================

int main(void) {
    printf("\n==========================================================\n");
    printf("   L-BRACKET METHOD OF MOMENTS SOLVER\n");
    printf("   Optimized SSOR with Hierarchical Omega Search\n");
    printf("==========================================================\n");

    clock_t program_start = clock();

    // Read configuration
    double strip_length, tx_x, tx_y;
    FILE *config = fopen("Transmitter pos", "r");
    if (!config || fscanf(config, "%lf %lf %lf", &strip_length, &tx_x, &tx_y) != 3) {
        fprintf(stderr, "ERROR: Cannot read 'Transmitter pos' file\n");
        fprintf(stderr, "Format: strip_length(m)  tx_x(m)  tx_y(m)\n");
        if (config) fclose(config);
        return 1;
    }
    fclose(config);

    Transmitter tx = {tx_x, tx_y};

    // Initialize problem
    EMParameters params = initialize_em_parameters(FREQUENCY_MHZ);
    Geometry geom = create_geometry(strip_length, params);

    printf("\nConfiguration:\n");
    printf("  Frequency: %.0f MHz (λ=%.4f m)\n", 
           params.frequency/1e6, params.wavelength);
    printf("  Strip length: %.4f m (%.2f λ)\n", 
           strip_length, strip_length/params.wavelength);
    printf("  Segments: N=%d (%d + %d)\n",
           geom.total_segments, geom.num_horizontal, geom.num_vertical);
    printf("  Transmitter: (%.2f, %.2f) m\n", tx.x, tx.y);
    printf("  OpenMP threads: %d\n", omp_get_max_threads());

    int N = geom.total_segments;
    double matrix_mb = (N * N * sizeof(double complex)) / (1024.0 * 1024.0);
    printf("  Memory: %.2f MB\n", matrix_mb);

    // Allocate system matrices
    double complex *Z = aligned_alloc(64, N * N * sizeof(double complex));
    double complex *V = aligned_alloc(64, N * sizeof(double complex));
    double complex *x = aligned_alloc(64, N * sizeof(double complex));

    if (!Z || !V || !x) {
        fprintf(stderr, "ERROR: Memory allocation failed\n");
        free_geometry(&geom);
        aligned_free(Z); aligned_free(V); aligned_free(x);
        return 1;
    }

    // STEP 1: Assemble system
    printf("\n=== STEP 1: Matrix Assembly ===\n");
    clock_t t1 = clock();
    assemble_system(Z, V, geom, tx, params);
    double time_assembly = (double)(clock() - t1) / CLOCKS_PER_SEC;
    printf("Assembly time: %.2f seconds\n", time_assembly);

    // STEP 2: Optimize omega
    printf("\n=== STEP 2: Omega Optimization ===\n");
    clock_t t2 = clock();
    double optimal_omega = optimize_omega_hierarchical(Z, V, N);
    double time_omega = (double)(clock() - t2) / CLOCKS_PER_SEC;
    printf("Optimization time: %.2f seconds\n", time_omega);

    // STEP 3: Solve system
    printf("\n=== STEP 3: Solve System ===\n");
    ConvergenceResult result = solve_ssor(Z, V, x, N, optimal_omega, 
                                         CONVERGENCE_TOL, MAX_ITERATIONS);

    double total_time = (double)(clock() - program_start) / CLOCKS_PER_SEC;

    // Print summary
    printf("\n==========================================================\n");
    printf("   PERFORMANCE SUMMARY\n");
    printf("==========================================================\n");
    printf("Matrix assembly:     %.2f seconds\n", time_assembly);
    printf("Omega optimization:  %.2f seconds\n", time_omega);
    printf("SSOR solver:         %.2f seconds (%d iterations)\n", 
           result.time_seconds, result.iterations);
    printf("----------------------------------------------------------\n");
    printf("Total time:          %.2f seconds\n", total_time);
    printf("==========================================================\n");

    // Save results
    save_results("current_distribution.txt", x, N, params, geom, tx,
                optimal_omega, result, total_time);
    save_geometry("geometry.txt", geom, params);

    printf("\nResults saved to:\n");
    printf("  - current_distribution.txt\n");
    printf("  - geometry.txt\n");

    // Cleanup
    free_geometry(&geom);
    aligned_free(Z); aligned_free(V); aligned_free(x);

    printf("\n*** Computation complete ***\n\n");
    return 0;
}