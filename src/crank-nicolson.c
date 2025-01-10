#include <stdio.h>
#include <complex.h> // library for complex numbers

double norm_squared(int N, double complex psi[N], double dx) {
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        res += psi[i]*conj(psi[i]);
    }
    // res *= dx;
    return res;
}

double potential(double x) {
    return 0.0;
}

void mul_tridiagmat_vec(int N, double complex A[N][N], double complex v[N]) {
    double complex res[N];

    // since the matrix is tridiagonal, the computation is more efficient -> O(N)
    res[0] = A[0][0] * v[0] + A[0][1] * v[1];
    for (int i = 1; i < N-1; i++) {
        res[i] = A[i][i-1] * v[i-1] + A[i][i] * v[i] + A[i][i+1] * v[i+1];
    }
    res[N-1] = A[N-1][N-2] * v[N-2] + A[N-1][N-1] * v[N-1];

    for (int i = 0; i < N; i++) {
        v[i] = res[i];
    }
    return;
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
    // fprintf(f_psi, "t,re,im,norm_sq\n");
    fprintf(f_psi, "t");
    for (int i = 0; i < N; i++) {
        fprintf(f_psi, ",re%d,im%d", i, i);
    }
    fprintf(f_psi, ",norm_sq\n");
    // fprintf(f_psi, "%.10f,%.10f,%.10f,%.10f\n", -dt, 1.0, 0.0, 1.0); // initial normalization

    // potential
    double V[N];
    for (int i = 0; i < N; i++) {
        V[i] = potential(-L + i * dx);
    }

    // explicit evolution matrix
    double complex A_half[N][N];
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

        double normalization = norm_squared(N, psi, dx);
        if (k % 10 == 0) {
            fprintf(f_psi, "%.10f", k*dt);
            for (int i = 0; i < N; i++) {
                fprintf(f_psi, ",%.10f,%.10f", creal(psi[i]), cimag(psi[i]));
            }
            fprintf(f_psi, ",%.10f\n", normalization);
        }
    }


    fclose(f_psi);
    return 0;
}