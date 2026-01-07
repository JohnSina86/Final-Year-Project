#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define MAX_ITER 50000
#define TOL 1e-9

// SSOR iteration function
int ssor_solve(double complex **Z, double complex *V, double complex *x, 
               int N, double omega, double tol, int max_iter) {

    double complex *x_new = malloc(N * sizeof(double complex));
    if (!x_new) {
        fprintf(stderr, "Memory allocation failed for x_new\n");
        return -1;
    }

    // Initial guess: diagonal preconditioning
    for (int i = 0; i < N; i++) {
        x[i] = V[i] / Z[i][i];
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

        for (int i = 0; i < N; i++) {
            x[i] = x_new[i];
        }

        if (relative_error < tol) {
            converged_iter = iter + 1;
            printf("Converged: %d iterations, error=%.2e\n", converged_iter, relative_error);
            break;
        }

        // Progress reporting every 1000 iterations
        if ((iter + 1) % 1000 == 0) {
            printf("  Iteration %d: error=%.2e\n", iter + 1, relative_error);
        }
    }

    free(x_new);
    return converged_iter;
}

// Simplified omega finder - single coarse pass
double find_omega(double complex **Z, double complex *V, int N) {
    printf("\nOptimizing relaxation parameter omega...\n");

    double best_omega = 1.5;
    double best_error = 1e20;
    int best_iters = MAX_ITER;

    // Quick coarse scan with loose tolerance
    for (double omega = 1.0; omega <= 1.95; omega += 0.05) {
        double complex *x_test = malloc(N * sizeof(double complex));
        double complex *x_new = malloc(N * sizeof(double complex));

        // Initial guess
        for (int i = 0; i < N; i++) {
            x_test[i] = V[i] / Z[i][i];
            x_new[i] = x_test[i];
        }

        // Quick test: 1000 iterations max, loose tolerance
        double final_error = 1e20;
        for (int iter = 0; iter < 1000; iter++) {
            // Forward sweep
            for (int i = 0; i < N; i++) {
                double complex sum = V[i];
                for (int j = 0; j < i; j++) sum -= Z[i][j] * x_new[j];
                sum -= Z[i][i] * x_test[i];
                for (int j = i + 1; j < N; j++) sum -= Z[i][j] * x_test[j];
                x_new[i] = x_test[i] + omega * (sum / Z[i][i]);
            }

            // Backward sweep
            for (int i = N - 1; i >= 0; i--) {
                double complex sum = V[i];
                for (int j = 0; j < i; j++) sum -= Z[i][j] * x_new[j];
                sum -= Z[i][i] * x_new[i];
                for (int j = i + 1; j < N; j++) sum -= Z[i][j] * x_new[j];
                x_new[i] = x_new[i] + omega * (sum / Z[i][i]);
            }

            // Check convergence
            double maxdiff = 0.0, maxval = 1e-20;
            for (int i = 0; i < N; i++) {
                maxdiff = fmax(maxdiff, cabs(x_new[i] - x_test[i]));
                maxval = fmax(maxval, cabs(x_new[i]));
            }
            final_error = maxdiff / maxval;

            for (int i = 0; i < N; i++) x_test[i] = x_new[i];

            if (final_error < 1e-3) break;
        }

        if (final_error < best_error) {
            best_error = final_error;
            best_omega = omega;
        }

        free(x_test);
        free(x_new);
    }

    printf("Optimal omega: %.4f (test error: %.2e)\n", best_omega, best_error);
    return best_omega;
}

void solve_mom_problem(double length_of_strip, double TX, double TY) {
    /* Physical constants */
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

    /* Geometry */
    double strip_length = length_of_strip;
    double angle = 0.0;
    double strip_start_x = 0.0 * lambda;
    double strip_start_y = 0.0;

    /* Discretization */
    int disc_per_lambda = 10;
    double delta_s = lambda / disc_per_lambda;

    /* Allocate position arrays */
    int max_N = 100000;
    double *x_pos = malloc(max_N * sizeof(double));
    double *y_pos = malloc(max_N * sizeof(double));

    if (!x_pos || !y_pos) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(x_pos); free(y_pos);
        return;
    }

    /* Calculate basis function positions */
    int N_per_strip = (int)ceil(strip_length / delta_s);
    if (N_per_strip > max_N) {
        fprintf(stderr, "Error: N_per_strip=%d exceeds max_N=%d\n", N_per_strip, max_N);
        free(x_pos); free(y_pos);
        return;
    }

    for (int ct = 0; ct < N_per_strip; ct++) {
        x_pos[ct] = strip_start_x + (ct + 0.5) * cos(angle) * delta_s;
        y_pos[ct] = strip_start_y + (ct + 0.5) * sin(angle) * delta_s;
    }

    int N = N_per_strip;
    printf("\n=== MoM Problem Setup ===\n");
    printf("Problem size: N = %d\n", N);
    printf("Strip length: %.4f lambda\n", strip_length / lambda);
    printf("Transmitter position: (%.4f, %.4f)\n", TX, TY);

    /* Allocate Z and V */
    double complex **Z = malloc(N * sizeof(double complex *));
    for (int i = 0; i < N; i++) {
        Z[i] = malloc(N * sizeof(double complex));
    }
    double complex *V = malloc(N * sizeof(double complex));
    double complex *x = malloc(N * sizeof(double complex));

    if (!V || !Z || !x) {
        fprintf(stderr, "Memory allocation failed for Z, V, or x.\n");
        free(x_pos); free(y_pos); free(V); free(x);
        for (int i = 0; i < N; i++) free(Z[i]);
        free(Z);
        return;
    }

    printf("\nBuilding impedance matrix Z and excitation vector V...\n");

    /* Build incident field vector V */
    for (int ct = 0; ct < N; ct++) {
        double Rx = x_pos[ct] - TX;
        double Ry = y_pos[ct] - TY;
        double R = sqrt(Rx * Rx + Ry * Ry);
        double J0 = j0(k0 * R);
        double Y0 = y0(k0 * R);
        V[ct] = J0 - I * Y0;
    }

    /* Diagonal self term */
    double complex self = -1.0 * delta_s *
                          (1.0 - I * (2.0 / MPI) *
                           log(1.781 * k0 * delta_s / (4.0 * exp(1.0))));

    /* Build impedance matrix Z */
    for (int ct1 = 0; ct1 < N; ct1++) {
        for (int ct2 = 0; ct2 < N; ct2++) {
            double complex val;
            if (ct1 != ct2) {
                double Rx = x_pos[ct1] - x_pos[ct2];
                double Ry = y_pos[ct1] - y_pos[ct2];
                double R = sqrt(Rx * Rx + Ry * Ry);
                double J0 = j0(k0 * R);
                double Y0 = y0(k0 * R);
                double complex the_hank = -(J0 - I * Y0);
                val = (k0 * eta / 4.0) * the_hank * delta_s;
            } else {
                val = (k0 * eta / 4.0) * self;
            }
            Z[ct1][ct2] = val;
        }
    }

    printf("Matrix assembly complete.\n");

    /* Find optimal omega (comment out to use fixed omega=1.7) */
    double optimal_omega = find_omega(Z, V, N);
    // double optimal_omega = 1.7;  // Uncomment to skip optimization

    /* Solve using SSOR */
    printf("\n=== Solving with SSOR (omega=%.4f) ===\n", optimal_omega);
    int iters = ssor_solve(Z, V, x, N, optimal_omega, TOL, MAX_ITER);

    if (iters == -1) {
        printf("Warning: SSOR did not converge within %d iterations\n", MAX_ITER);
    } else {
        printf("\nSolution successful!\n");
    }

    /* Save result */
    FILE *fx = fopen("current_distribution.txt", "w");
    if (fx) {
        fprintf(fx, "%% Current distribution solution (SSOR method)\n");
        fprintf(fx, "%% N=%d, omega=%.4f, iterations=%d\n", N, optimal_omega, iters);
        for (int i = 0; i < N; i++) {
            fprintf(fx, "(%15.8e,%15.8e)\n", creal(x[i]), cimag(x[i]));
        }
        fclose(fx);
        printf("Current distribution saved to current_distribution.txt\n");
    }

    /* Cleanup */
    free(x_pos);
    free(y_pos);
    free(V);
    free(x);
    for (int i = 0; i < N; i++) free(Z[i]);
    free(Z);
}

int main(void) {
    double length_of_strip, TX, TY;

    FILE *f_pos = fopen("Transmitter pos", "r");
    if (f_pos == NULL) {
        fprintf(stderr, "Error opening 'Transmitter pos' file.\n");
        return 1;
    }
    fscanf(f_pos, "%lf", &length_of_strip);
    fscanf(f_pos, "%lf", &TX);
    fscanf(f_pos, "%lf", &TY);
    fclose(f_pos);

    solve_mom_problem(length_of_strip, TX, TY);
    return 0;
}
