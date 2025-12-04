#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void solve_complex_system(double complex *Z, double complex *V, double complex *x, int N) {
    // Create augmented matrix [Z | V]
    double complex **aug = malloc(N * sizeof(double complex *));
    for (int i = 0; i < N; i++) {
        aug[i] = malloc((N + 1) * sizeof(double complex));
        for (int j = 0; j < N; j++) {
            aug[i][j] = Z[i * N + j];
        }
        aug[i][N] = V[i];
    }

    // Gaussian elimination with partial pivoting
    for (int k = 0; k < N; k++) {
        // Find pivot
        int pivot = k;
        double max_val = cabs(aug[k][k]);
        for (int i = k + 1; i < N; i++) {
            if (cabs(aug[i][k]) > max_val) {
                max_val = cabs(aug[i][k]);
                pivot = i;
            }
        }

        // Swap rows if needed
        if (pivot != k) {
            double complex *temp = aug[k];
            aug[k] = aug[pivot];
            aug[pivot] = temp;
        }

        // Forward elimination
        for (int i = k + 1; i < N; i++) {
            double complex factor = aug[i][k] / aug[k][k];
            for (int j = k; j <= N; j++) {
                aug[i][j] -= factor * aug[k][j];
            }
        }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
        x[i] = aug[i][N];
        for (int j = i + 1; j < N; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }

    // Free augmented matrix
    for (int i = 0; i < N; i++) {
        free(aug[i]);
    }
    free(aug);
}

void make_the_mom_matrix(double length_of_strip, double TX, double TY) {
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
    int no_of_strips = 1;
    double strip_length = length_of_strip;
    double angle = 0.0;

    double strip_start_x = 0.0 * lambda;
    double strip_start_y = 0.0;

    /* Discretisation */
    int disc_per_lambda = 10;
    double delta_s = lambda / disc_per_lambda;

    /* For simplicity, assume a max number of basis functions */
    int max_N = 10000;
    double *x_pos = (double *)malloc(max_N * sizeof(double));
    double *y_pos = (double *)malloc(max_N * sizeof(double));
    double *arc_length = (double *)malloc(max_N * sizeof(double));

    if (!x_pos || !y_pos || !arc_length) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(x_pos); free(y_pos); free(arc_length);
        return;
    }

    /* Count basis functions */
    int pos_counter = 0;
    int N_per_strip = (int)ceil(strip_length / delta_s);

    if (N_per_strip > max_N) {
        fprintf(stderr, "Increase max_N, N_per_strip=%d\n", N_per_strip);
        free(x_pos); free(y_pos); free(arc_length);
        return;
    }

    for (int ct_within_strip = 1; ct_within_strip <= N_per_strip; ++ct_within_strip) {
        pos_counter++;
        int idx = pos_counter - 1;
        x_pos[idx] = strip_start_x + (ct_within_strip - 0.5) * cos(angle) * delta_s;
        y_pos[idx] = strip_start_y + (ct_within_strip - 0.5) * sin(angle) * delta_s;
        arc_length[idx] = (pos_counter - 0.5) * delta_s;
    }

    int N = pos_counter;
    printf("Problem size: N = %d\n", N);

    /* Allocate V and Z */
    double complex *V = (double complex *)malloc(N * sizeof(double complex));
    double complex *Z = (double complex *)malloc(N * N * sizeof(double complex));
    if (!V || !Z) {
        fprintf(stderr, "Memory allocation failed for V or Z.\n");
        free(x_pos); free(y_pos); free(arc_length); free(V); free(Z);
        return;
    }

    /* Incident field */
    for (int ct = 0; ct < N; ++ct) {
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

    /* Fill Z */
    for (int ct1 = 0; ct1 < N; ++ct1) {
        for (int ct2 = 0; ct2 < N; ++ct2) {
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
            Z[ct1 * N + ct2] = val;
        }
    }

    /* Solve Z*x = V for exact solution */
    printf("Solving for exact current distribution x...\n");
    double complex *x_exact = (double complex *)malloc(N * sizeof(double complex));
    if (!x_exact) {
        fprintf(stderr, "Memory allocation failed for x_exact.\n");
        free(x_pos); free(y_pos); free(arc_length); free(V); free(Z);
        return;
    }

    solve_complex_system(Z, V, x_exact, N);
    printf("Solution complete.\n");

    /* Write exact solution to file */
    FILE *fx = fopen("x_exact.txt", "w");
    if (fx) {
        fprintf(fx, "%% Exact solution x (current coefficients)\n");
        for (int i = 0; i < N; i++) {
            fprintf(fx, "(%15.8e,%15.8e)\n", creal(x_exact[i]), cimag(x_exact[i]));
        }
        fclose(fx);
        printf("Exact solution written to x_exact.txt\n");
    }

    /* Also write to old format for compatibility */
    FILE *fp = fopen("function_result.txt", "w");
    if (!fp) {
        fprintf(stderr, "Could not open output file.\n");
        free(x_pos); free(y_pos); free(arc_length); free(V); free(Z); free(x_exact);
        return;
    }

    fprintf(fp, "%% Z matrix (rows = basis functions, columns = basis functions)\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double complex zval = Z[i * N + j];
            fprintf(fp, "(%15.8e,%15.8e) ", creal(zval), cimag(zval));
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n%% V vector (incident field)\n");
    for (int i = 0; i < N; i++) {
        fprintf(fp, "(%15.8e,%15.8e)\n", creal(V[i]), cimag(V[i]));
    }

    fclose(fp);

    free(x_pos);
    free(y_pos);
    free(arc_length);
    free(V);
    free(Z);
    free(x_exact);
}

int main(void) {
    // Try these parameters to get a different current distribution
    double length_of_strip;
    double TX;
    double TY;

    FILE *f_pos = fopen("Transmitter pos", "r");
    if (f_pos == NULL) {
        fprintf(stderr, "Error opening 'Transmitter pos' file.\n");
        return 1;
    }
    fscanf(f_pos, "%lf", &length_of_strip);
    fscanf(f_pos, "%lf", &TX);
    fscanf(f_pos, "%lf", &TY);
    fclose(f_pos);
    
    make_the_mom_matrix(length_of_strip, TX, TY);
    return 0;
}

