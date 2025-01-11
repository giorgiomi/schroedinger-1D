#include <stdio.h>
#include <complex.h> // library for complex numbers
#include <stdlib.h>
#include "functions.h"

double potential(double x) {
    return 0.0;
}

int main(int arcv, char** argv) {
    // parameters
    double L = 1.0;                             // box size
    int N = 100;                                // number of grid separations
    int M = 30000;                              // number of time steps
    double dx = 2 * L / (double)(N + 1);        // space interval
    double dt = 1e-6;                           // time interval
    double complex dtau = - dt * I;             // complex tau interval
    double complex eta = - dtau / (dx * dx);    // eta parameter

    // files
    FILE* f_psi = fopen("data/C-N.csv", "w");
    fprintf(f_psi, "t");
    for (int i = 0; i < N; i++) {
        fprintf(f_psi, ",re%d,im%d", i, i);
    }
    fprintf(f_psi, ",norm_sq\n");

    // potential
    double V[N];
    for (int i = 0; i < N; i++) {
        V[i] = potential(-L + (i + 1) * dx);
    }

    // explicit evolution matrix (half time step)
    double complex **A_half = malloc(N * sizeof(_Complex double *));
    for (int i = 0; i < N; i++) {
        A_half[i] = malloc(N * sizeof(_Complex double));
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                A_half[i][j] = 1.0 - eta + dtau * V[i] / 2.0;
            } else if (i == j + 1 || j == i + 1) {
                A_half[i][j] = eta / 2.0;
            } else {
                A_half[i][j] = 0.0;
            }
        }
    }

    // implicit C-N evolution coefficients
    double complex a[N];
    double complex alpha[N];
    double complex beta[N];
    double complex gamma = -eta/2;
    
    a[0] = 1 + eta - V[0] * dtau / 2;
    alpha[0] = a[0];
    for (int i = 1; i < N; i++) {
        a[i] = 1 + eta - V[i] * dtau / 2;
        beta[i] = -eta/(2*alpha[i-1]);
        alpha[i] = a[i] - beta[i]/2;
    }
    

    // initial condition
    double complex psi[N];
    psi[0] = 1.0;
    for (int i = 1; i < N; i++) {
        psi[i] = 0.0;
    }

    // simulation
    double complex y[N];
    for (int k = 0; k < M; k++) {
        // explicit half-step
        mul_tridiagmat_vec(N, A_half, psi);

        // implicit half-step
        y[0] = psi[0];
        for (int i = 1; i < N; i++) {
            y[i] = psi[i] - beta[i] * y[i-1];
        }
        psi[N-1] = y[N-1]/alpha[N-1];
        for (int i = N-2; i >= 0; i--) {
            psi[i] = y[i] / alpha[i] - psi[i+1] * gamma / alpha[i];
        }
        // evolution done!

        
        if (k % 10 == 0) {
            double normalization = norm_squared(N, psi, dx);
            printLineOnFile(f_psi, N, k*dt, psi, normalization);
        }
    }

    for (int i = 0; i < N; i++) {
        free(A_half[i]);
    }
    free(A_half);
    fclose(f_psi);
    return 0;
}