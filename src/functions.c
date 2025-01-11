#include <stdio.h>
#include <complex.h>

double normSquared(int N, double complex *psi, double dx) {
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        res += psi[i]*conj(psi[i]);
    }
    // res *= dx;
    return res;
}

double normSquaredSimpson(int N, double complex *psi, double dx){
    double integral = (psi[N-1] * conj(psi[N-1])) * dx / 3;

    // assume N even
    for (int i = 1; i <= N/2 - 1; i++) {
        integral += (2 * dx / 3) * (psi[2 * i - 1] * conj(psi[2 * i - 1]));
        integral += (4 * dx / 3) * (psi[2 * i - 2] * conj(psi[2 * i - 2]));
    }
    integral += (4 * dx / 3) * (psi[N-2] * conj(psi[N-2]));

    return integral;
}

void mul_tridiagmat_vec(int N, double complex **A, double complex *v) {
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

void printLineOnFile(FILE *f, int N, double time, double _Complex *psi, double norm) {
    fprintf(f, "%.10f", time);
        for (int i = 0; i < N; i++) {
            fprintf(f, ",%.10f,%.10f", creal(psi[i]), cimag(psi[i]));
        }
    fprintf(f, ",%.10f\n", norm);
    return;
}