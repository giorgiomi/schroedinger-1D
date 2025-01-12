#include <stdio.h>
#include <complex.h>

double normSquared(int N, double *f, double dx) {
    double res = 0.0;
    for (int i = 0; i < N; i++) {
        res += f[i];
    }
    res *= dx;
    return res;
}

double normSquaredLeft(int N, double *f, double dx) {
    double res = 0.0;
    for (int i = 0; i <= (N - 1)/2; i++) {
        res += f[i];
    }
    res *= dx;
    return res;
}

double normSquaredRight(int N, double *f, double dx) {
    double res = 0.0;
    for (int i = (N + 1)/2; i < N; i++) {
        res += f[i];
    }
    res *= dx;
    return res;
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

void printLineOnFile(FILE *f, int N, double time, double _Complex *psi, double norm, double x_mean, double x2_mean) {
    fprintf(f, "%.10f", time);
        for (int i = 0; i < N; i++) {
            fprintf(f, ",%.10f,%.10f", creal(psi[i]), cimag(psi[i]));
        }
    fprintf(f, ",%.10f,%.10f,%.10f\n", norm, x_mean, x2_mean);
    return;
}