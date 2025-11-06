// mom_simple.c
// Minimal 2D MoM matrix builder for one straight strip.
// - Discretize strip into N segments
// - Z_ij ~ i*omega*mu0 * (Δs/4) * H0^(2)(k0 * R_ij)
// - V_i  ~ H0^(2)(k0 * R_src_i)
// Uses j0() and y0() to make Hankel: H0^(2) = J0 - i Y0.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Return H0^(2)(x) = J0(x) - i*Y0(x)
static inline double complex hankel0_2(double x) {
    // j0, y0 are POSIX (math.h). If unavailable on your platform,
    // replace with approximations or link to a special function lib.
    return j0(x) - I * y0(x);
}

int main(void) {
    // --------- Inputs (tweak as you like) ----------
    double length_of_strip = 1.0;   // meters
    double TX = 0.30, TY = 0.40;    // point source location (m)

    // --------- Physical constants ----------
    const double f = 1.0e9;                   // 1 GHz
    const double epsilon0 = 8.854187817e-12;  // F/m
    const double mu0 = 4.0e-7 * M_PI;         // H/m
    const double c0 = 1.0 / sqrt(epsilon0 * mu0);
    const double omega = 2.0 * M_PI * f;
    const double k0 = omega / c0;             // wavenumber
    (void)c0; // silence unused warning if not printed

    // --------- Discretization ----------
    const int   disc_per_lambda = 20;               // segments per wavelength
    const double lambda = 2.0 * M_PI / k0;
    int N = (int)ceil((length_of_strip / lambda) * disc_per_lambda);
    if (N < 1) N = 1;

    const double delta_s = length_of_strip / (double)N; // uniform segment length

    // Segment centers (straight strip along +x, starting at x=0, y=0)
    double *x = (double*)malloc(N * sizeof(double));
    double *y = (double*)malloc(N * sizeof(double));
    if (!x || !y) { fprintf(stderr, "alloc failed\n"); return 1; }

    for (int i = 0; i < N; ++i) {
        x[i] = (i + 0.5) * delta_s; // center of segment i
        y[i] = 0.0;
    }

    // --------- Build Z (NxN complex) and V (N complex) ----------
    double complex *Z = (double complex*)calloc((size_t)N * N, sizeof(double complex));
    double complex *V = (double complex*)calloc((size_t)N, sizeof(double complex));
    if (!Z || !V) { fprintf(stderr, "alloc failed\n"); return 1; }

    // RHS: incident field from point source (simple model)
    // V_i ~ H0^(2)(k0 * R_src_i). (Scaling constants omitted for simplicity.)
    for (int i = 0; i < N; ++i) {
        double dx = x[i] - TX;
        double dy = y[i] - TY;
        double Rsrc = hypot(dx, dy);
        if (Rsrc < 1e-12) Rsrc = 1e-12; // avoid singularity
        V[i] = hankel0_2(k0 * Rsrc);
    }

    // Impedance: off-diagonal uses true distance; diagonal uses R ~ Δs/2
    // Z_ij ≈ i*ω*μ0 * (Δs/4) * H0^(2)(k0 * R_ij)
    const double complex pref = I * omega * mu0 * (delta_s / 4.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double Rij;
            if (i == j) {
                // crude self-term regularization
                Rij = delta_s * 0.5;
            } else {
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                Rij = hypot(dx, dy);
            }
            if (Rij < 1e-12) Rij = 1e-12;
            Z[(size_t)i * N + j] = pref * hankel0_2(k0 * Rij);
        }
    }

    // --------- Output ----------
    printf("N = %d, delta_s = %.6g m, lambda = %.6g m\n", N, delta_s, lambda);
    printf("\n-- V (RHS) --\n");
    for (int i = 0; i < N; ++i) {
        printf("V[%d] = %.6e%+.6ei\n", i, creal(V[i]), cimag(V[i]));
    }

    printf("\n-- Z (Impedance Matrix) --\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double complex zij = Z[(size_t)i * N + j];
            printf("(% .3e%+.3ei)%s", creal(zij), cimag(zij), (j==N-1) ? "" : " ");
        }
        printf("\n");
    }

    // Clean up
    free(x); free(y); free(Z); free(V);
    return 0;
}
