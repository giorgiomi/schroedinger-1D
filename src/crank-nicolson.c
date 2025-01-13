#include <stdio.h>
#include <complex.h> // library for complex numbers
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double potential(double x, double V0, double a) {
    return V0 * (x * x - a) * (x * x - a);
    // return 0.0;
}

int main(int argc, char** argv) {
    // parameters
    int N = 199;                                    // number of grid separations
    int M = 4.0e3;                                  // number of time steps
    double L = 1.0;                                 // box size
    double dx = 2 * L / (double)(N + 1);            // space interval
    double dt = 1e-4;                               // time interval
    double complex dtau = - dt * I;                 // complex tau interval
    double complex eta = - dtau / (2 * dx * dx);    // eta parameter
    double V0;                                      // potential strength
    double A;                                       // double well separation parameter
    int n_print = 1;                               // print on file every n_print iteration

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <V0> <A>\n", argv[0]);
        return 1;
    }
    V0 = atof(argv[1]);
    A = atof(argv[2]);

    // printf("==================================================================================\n");    
    // printf("Running C-N with N = %d, M = %d, L = %.2f, dx = %.4e, dt = %.2e\n\n", N, M, L, dx, dt);

    // files
    FILE* f_psi = fopen("data/trapped/C-N.csv", "w");
    printHeaderOnFile(f_psi, N, 0);

    FILE* f_param = fopen("data/trapped/param.csv", "w");
    fprintf(f_param, "N,M,L,dx,dt,V0,a\n");
    fprintf(f_param, "%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f\n", N, M, L, dx, dt, V0, A);

    // potential
    double V[N];
    for (int i = 0; i < N; i++) {
        V[i] = potential(-L + (i + 1) * dx, V0, A);
        // printf("V[%d] = %.4f\n", i, V[i]);
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
        beta[i] = -eta / (2 * alpha[i-1]);
        alpha[i] = a[i] - (eta * eta) / (4 * alpha[i-1]);
    }
    
    // initial condition
    double complex psi[N];
    int i_start = (int)((L - sqrt(A)) / dx) - 1;
    // printf("i_start = %d", i_start);
    psi[i_start] = 1.0 / sqrt(dx);
    for (int i = 0; i < N; i++) {
        if (i != i_start) {
            psi[i] = 0.0;
        }
    }

    // simulation
    for (int k = 0; k < M; k++) {
        // printf("\rStep %d of %d", k + 1, M);
        // fflush(stdout);

        // evolution step in one line
        crankNicolsonStep(N, A_half, psi, alpha, beta, gamma);
        
        // printing on file every n_print iterations 
        if (k % n_print == 0) {
            double psi_norm[N];
            for (int i = 0; i < N; i++) {
                psi_norm[i] = psi[i] * conj(psi[i]);
            }
            double prob_left = normSquaredLeft(N, psi_norm, dx);
            double prob_right = normSquaredRight(N, psi_norm, dx);
            printLineOnFile(f_psi, N, k*dt, NULL, 1.0, prob_left, prob_right);
        }
    }

    // printf("\n\nSimulation completed 🎉🎊\n");
    // printf("==================================================================================\n");

    for (int i = 0; i < N; i++) {
        free(A_half[i]);
    }
    free(A_half);
    fclose(f_psi);
    return 0;
}