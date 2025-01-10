#include <stdio.h>
#include <complex.h> // library for complex numbers

double norm_squared(int N, double complex psi[N], double dx) {
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        res += psi[i]*conj(psi[i]);
    }
    res += dx;
    return res;
}

double V(double x) {
    return 0.0;
}

void mul_mat_vec(int N, double complex A[N][N], double complex v[N]) {
    double complex res[N];
    for (int i = 0; i < N; i++) {
        res[i] = 0.0;
        for (int j = 0; j < N; j++) {
            res[i] += A[i][j] * v[j]; //this can be made more efficient, one needs to calculate only three values!
        }
    }
    for (int i = 0; i < N; i++) {
        v[i] = res[i];
    }
    return;
}

int main(int arcv, char** argv) {
    // parameters
    double L = 1.0;                         // box size
    int N = 100;                            // number of grid separations
    int M = 10000;                           // number of time steps
    double dx = 2 * L / (double)(N + 1);    // space interval
    double dt = 1e-6;                        // time interval
    double complex dtau = - dt * I;         // complex tau interval
    double complex eta = - dtau / (dx * dx);

    // files
    FILE* f_psi = fopen("data/psi.csv", "w");
    fprintf(f_psi, "t,re,im,norm_sq\n");
    // fprintf(f_psi, "%.10f,%.10f,%.10f,%.10f\n", -dt, 1.0, 0.0, 1.0); // initial normalization

    // evolution matrix
    double complex A[N][N];
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
        // printf("%d\n", k);
        // printf("%.3f\n", creal(psi[0]));
        mul_mat_vec(N, A, psi);             // evolution step
        double normalization = norm_squared(N, psi, dx);
        fprintf(f_psi, "%.10f,%.10f,%.10f,%.10f\n", k*dt, creal(psi[50]), cimag(psi[50]), normalization);
    }


    fclose(f_psi);
    return 0;
}