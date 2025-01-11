#include <stdio.h>
#include <complex.h> // library for complex numbers
#include <stdlib.h>
#include "functions.h"

int main(int arcv, char** argv) {
    // parameters
    double L = 1.0;                             // box size
    int N = 50;                                // number of grid separations
    int M = 30000;                              // number of time steps
    double dx = 2 * L / (double)(N + 1);        // space interval
    double dt = 1e-6;                           // time interval
    double complex dtau = - dt * I;             // complex tau interval
    double complex eta = - dtau / (dx * dx);    // eta parameter

    // files
    FILE* f_psi = fopen("data/EU.csv", "w");
    fprintf(f_psi, "t");
    for (int i = 0; i < N; i++) {
        fprintf(f_psi, ",re%d,im%d", i, i);
    }
    fprintf(f_psi, ",norm_sq\n");

    // evolution matrix
    double complex **A = malloc(N * sizeof(_Complex double *));
    for (int i = 0; i < N; i++) {
        A[i] = malloc(N * sizeof(_Complex double));
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A[i][j] = 1.0 - 2.0 * eta;
            } else if (i == j + 1 || j == i + 1) {
                A[i][j] = eta;
            } else {
                A[i][j] = 0.0;
            }
        }
    }

    // initial condition
    double complex psi[N];
    psi[0] = 1.0;
    for (int i = 1; i < N; i++) {
        psi[i] = 0.0;
    }

    // simulation
    for (int k = 0; k < M; k++) {
        mul_tridiagmat_vec(N, A, psi);             // evolution step

        if (k % 10 == 0) {
            double normalization = normSquaredSimpson(N, psi, dx);
            printLineOnFile(f_psi, N, k*dt, psi, normalization);
        }
    }

    for (int i = 0; i < N; i++) {
        free(A[i]);
    }
    free(A);
    fclose(f_psi);
    return 0;
}