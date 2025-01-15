// SOLVES THE 1-D SCHROEDINGER EQUATION FOR A FINITE BOX AND A DOUBLE WELL POTENTIAL
// TROTTER-SUZUKI METHOD
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
    int M = 5.0e3;                                  // number of time steps
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

    printf("==================================================================================\n");    
    printf("Running T-S with N = %d, M = %d, L = %.2f, dx = %.4e, dt = %.2e\n\n", N, M, L, dx, dt);

    // files
    FILE* f_psi = fopen("data/test/T-S.csv", "w");
    printHeaderOnFile(f_psi, N, 1);

    FILE* f_energy = fopen("data/test/energy.csv", "w");
    fprintf(f_energy, "t,K,V,E\n");

    FILE* f_param = fopen("data/test/paramT-S.csv", "w");
    fprintf(f_param, "N,M,L,dx,dt,V0,a\n");
    fprintf(f_param, "%d,%d,%.10f,%.10f,%.10f,%.10f,%.10f\n", N, M, L, dx, dt, V0, A);

    // potential
    double V[N];
    for (int i = 0; i < N; i++) {
        V[i] = potential(-L + (i + 1) * dx, V0, A);
        // printf("V[%d] = %.4f\n", i, V[i]);
    }
    
    // initial condition gaussian wave packet
    double complex psi[N];
    double sigma = 1.0 * sqrt(dx);
    // double sigma = 0.001 / dx;
    for (int i = 0; i < N; i++) {
        double x = - L + (i + 1) * dx;
        psi[i] = wf(x, -sqrt(A),  pow(1/(2*M_PI*sigma*sigma), 0.25), sigma, 0.0);
    }
    printLineOnFile(f_psi, N, 0.0, psi, 1.0, 1.0, 0.0);

    // simulation
    for (int k = 0; k < M; k++) {
        printf("\rStep %d of %d", k + 1, M);
        fflush(stdout);

        // evolution step in one line
        trotterSuzukiStep(N, V, dt, dx, psi, eta);
        
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

    printf("\n\nSimulation completed 🎉🎊\n");
    printf("==================================================================================\n");
    fclose(f_psi);
    return 0;
}