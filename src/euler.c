#include <stdio.h>
#include <complex.h> // library for complex numbers
#include <stdlib.h>
#include <math.h>
#include "functions.h"

int main(int arcv, char** argv) {
    // parameters
    double L = 1.0;                             // box size
    int N = 199;                                // number of grid separations
    int M = 3e4;                              // number of time steps
    double dx = 2 * L / (double)(N + 1);        // space interval
    double dt = 1e-6;                           // time interval
    double complex dtau = - dt * I;             // complex tau interval
    double complex eta = - dtau / (dx * dx);    // eta parameter

    printf("==================================================================================\n");    
    printf("Running Eu with N = %d, M = %d, L = %.2f, dx = %.4e, dt = %.2e\n\n", N, M, L, dx, dt);

    // files
    FILE* f_psi = fopen("data/free/EU.csv", "w");
    fprintf(f_psi, "t");
    for (int i = 0; i < N; i++) {
        fprintf(f_psi, ",re%d,im%d", i, i);
    }
    fprintf(f_psi, ",norm_sq,x,x2\n");

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
    psi[(N - 1)/2] = 1.0;
    for (int i = 0; i < N; i++) {
        if (i != (N - 1) / 2) {
            psi[i] = 0.0;
        }
    }
    // for (int i = 0; i < N; i++) {
    //     psi[i] = cos(M_PI * (-L + i * dx) / 2);
    // }

    // simulation
    for (int k = 0; k < M; k++) {
        mul_tridiagmat_vec(N, A, psi);             // evolution step

        if (k % 10 == 0) {
            double psi_norm[N];
            double x_psi_norm[N];
            double x2_psi_norm[N];
            for (int i = 0; i < N; i++) {
                double x = -L + (i + 1) * dx;
                psi_norm[i] = psi[i] * conj(psi[i]);
                x_psi_norm[i] = x * psi_norm[i];
                x2_psi_norm[i] = x * x * psi_norm[i];
            }
            double normalization = normSquared(N, psi_norm, dx);
            double x_mean = normSquared(N, x_psi_norm, dx);
            double x2_mean = normSquared(N, x2_psi_norm, dx);
            printLineOnFile(f_psi, N, k*dt, psi, normalization, x_mean, x2_mean);
        }
    }

    printf("\n\nSimulation completed 🎉🎊\n");
    printf("==================================================================================\n");

    for (int i = 0; i < N; i++) {
        free(A[i]);
    }
    free(A);
    fclose(f_psi);
    return 0;
}