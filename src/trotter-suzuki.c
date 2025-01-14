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
    int N = 49;                                    // number of grid separations
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
    // printf("Running T-S with N = %d, M = %d, L = %.2f, dx = %.4e, dt = %.2e\n\n", N, M, L, dx, dt);

    // files
    FILE* f_psi = fopen("data/test/T-S.csv", "w");
    printHeaderOnFile(f_psi, N, 1);

    FILE* f_param = fopen("data/test/paramT-S.csv", "w");
    fprintf(f_param, "N,M,L,dx,dt,V0,a\n");
    fprintf(f_param, "%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f\n", N, M, L, dx, dt, V0, A);

    // potential
    double V[N];
    for (int i = 0; i < N; i++) {
        V[i] = potential(-L + (i + 1) * dx, V0, A);
        // printf("V[%d] = %.4f\n", i, V[i]);
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
        trotterSuzukiStep(N, V, dt, dx, psi, eta);
        
        // printing on file every n_print iterations 
        if (k % n_print == 0) {
            double psi_norm[N];
            for (int i = 0; i < N; i++) {
                psi_norm[i] = psi[i] * conj(psi[i]);
            }
            double prob_left = normSquaredLeft(N, psi_norm, dx);
            double prob_right = normSquaredRight(N, psi_norm, dx);
            printLineOnFile(f_psi, N, k*dt, psi, 1.0, prob_left, prob_right);
        }
    }

    // printf("\n\nSimulation completed 🎉🎊\n");
    // printf("==================================================================================\n");

    
    fclose(f_psi);
    return 0;
}