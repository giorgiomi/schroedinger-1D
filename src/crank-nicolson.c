#include <stdio.h>
#include <complex.h> // library for complex numbers
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double potential(double x) {
    return 0.0;
}

int main(int arcv, char** argv) {
    // parameters
    int N = 199;                                // number of grid separations
    int M = 3e4;                                // number of time steps
    double L = 1.0;                             // box size
    double dx = 2 * L / (double)(N + 1);        // space interval
    double dt = 1e-6;                           // time interval
    double complex dtau = - dt * I;             // complex tau interval
    double complex eta = - dtau / (dx * dx);    // eta parameter

    printf("==================================================================================\n");    
    printf("Running C-N with N = %d, M = %d, L = %.2f, dx = %.4e, dt = %.2e\n\n", N, M, L, dx, dt);

    // files
    FILE* f_psi = fopen("data/free/C-N.csv", "w");
    fprintf(f_psi, "t");
    for (int i = 0; i < N; i++) {
        fprintf(f_psi, ",re%d,im%d", i, i);
    }
    fprintf(f_psi, ",norm_sq,x,x2\n");

    FILE* f_param = fopen("data/free/param.csv", "w");
    fprintf(f_param, "N,M,L,dx,dt\n");
    fprintf(f_param, "%d,%d,%.10f,%.10f,%.10f\n", N, M, L, dx, dt);

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
        alpha[i] = a[i] - (eta*eta)/(4*alpha[i-1]);
    }
    

    // initial condition
    double complex psi[N];
    psi[(N - 1)/2] = 1.0 / sqrt(dx);
    for (int i = 0; i < N; i++) {
        if (i != (N - 1) / 2) {
            psi[i] = 0.0;
        }
    }
    // for (int i = 0; i < N; i++) {
    //     // psi[i] = cos(M_PI * (-L + i * dx) / 2);
    //     psi[i] = sin(M_PI * (-L + i * dx)*2);
    // }

    // simulation
    double complex y[N];
    for (int k = 0; k < M; k++) {
        printf("\rStep %d of %d", k + 1, M);
        fflush(stdout);
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
        free(A_half[i]);
    }
    free(A_half);
    fclose(f_psi);
    return 0;
}