// SOLVES THE 1-D SCHROEDINGER EQUATION FOR A FINITE BOX AND A DOUBLE WELL POTENTIAL
#include <stdio.h>
#include <complex.h> // library for complex numbers
#include <stdlib.h>
#include <math.h>
#include "functions.h"

double potential(double x, double V0, double a) {
    return V0 * (x * x - a) * (x * x - a);
    // return 0.0;
}

double complex wf(double x, double x0, double A, double sigma, double k) {
    return A * exp(-(x - x0) * (x - x0)/(4*sigma*sigma))*cexp(I * k * x);
}

int main(int argc, char** argv) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <N> <dt> <V0> <A>\n", argv[0]);
        return 1;
    }

    // parameters
    int N = atoi(argv[1]);                          // number of grid separations
    int M = 1.5e3;                                  // number of time steps
    double L = 1.0;                                 // box size
    double dx = 2 * L / (double)(N + 1);            // space interval
    double dt = atof(argv[2]);                      // time interval
    // double dt = dx * dx;                            // time interval to keep same eta
    double complex dtau = - dt * I;                 // complex tau interval
    double complex eta = - dtau / (2 * dx * dx);    // eta parameter
    // printf("eta = %.2f + %.2fi\n", creal(eta), cimag(eta));
    double V0 = atof(argv[3]);                      // potential strength
    double A = atof(argv[4]);                       // double well separation parameter
    int n_print = 1;                                // print on file every n_print iteration

    // printf("==================================================================================\n");    
    // printf("Running C-N with N = %d, M = %d, L = %.2f, dx = %.4e, dt = %.2e\n\n", N, M, L, dx, dt);

    // files
    FILE* f_psi = fopen("data/trapped/C-N.csv", "w");
    printHeaderOnFile(f_psi, N, 1);

    FILE* f_energy = fopen("data/trapped/energy.csv", "w");
    fprintf(f_energy, "t,K,V,E\n");

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

    // DEBUG: PRINT A_HALF
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         printf("(%2.1f + %2.6fi) ", creal(A_half[i][j]), cimag(A_half[i][j]));
    //     }
    //     printf("\n");
    // }

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
    // int i_start = (int)((L - sqrt(A)) / dx) - 1;
    // psi[i_start] = 1.0 / sqrt(dx);
    // for (int i = 0; i < N; i++) {
    //     if (i != i_start) {
    //         psi[i] = 0.0;
    //     }
    // }

    // trying another initial condition
    double sigma = 1.0 * sqrt(dx);
    for (int i = 0; i < N; i++) {
        double x = - L + (i + 1) * dx;
        psi[i] = wf(x, -sqrt(A),  sqrt(1.0/(sigma * sqrt(2 * M_PI))), sigma, 0.0);
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
            double potential_energy[N];
            double kinetic_energy[N];
            kinetic_energy[0] = conj(psi[0]) * (- (psi[1] - 2 * psi[0]) / (2 * dx * dx));
            kinetic_energy[N-1] = conj(psi[N-1]) * (- (psi[N-2] - 2 * psi[0]) / (2 * dx * dx)); 
            for (int i = 0; i < N; i++) {
                psi_norm[i] = psi[i] * conj(psi[i]);
                potential_energy[i] = V[i] * psi_norm[i];
                if (i != 0 && i != N-1) {
                    kinetic_energy[i] = conj(psi[i]) * (- (psi[i+1] - 2 * psi[i] + psi[i-1]) / (2 * dx * dx));
                }
            }
            double V_avg = normSquared(N, potential_energy, dx);
            double K_avg = normSquared(N, kinetic_energy, dx);
            double prob_left = normSquaredLeft(N, psi_norm, dx);
            double prob_right = normSquaredRight(N, psi_norm, dx);
            printLineOnFile(f_psi, N, k*dt, psi, 1.0, prob_left, prob_right);
            printLineOnFile(f_energy, N, k*dt, NULL, K_avg, V_avg, K_avg + V_avg);
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