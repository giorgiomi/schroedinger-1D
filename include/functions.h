#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <complex.h>

/*
Computes the normalization of a function at a fixed time
N: array length
f: function array
dx: spatial step size
*/
double normSquared(int N, double *f, double dx);

/*
Computes the half-integral (L) of a function at a fixed time
N: array length
f: function array
dx: spatial step size
*/
double normSquaredLeft(int N, double *f, double dx);

/*
Computes the half-integral (R) of a function at a fixed time
N: array length
f: function array
dx: spatial step size
*/
double normSquaredRight(int N, double *f, double dx);



/*
Computes the multiplication between a tridiagonal matrix and a vector
N: vector length
A: tridiagonal matrix
v: vector
*/
void mul_tridiagmat_vec(int N, double _Complex **A, double _Complex *v);

/*
Prints line on file
f: file pointer
N: wave function array length
time: specify time to print
psi: wave function array
norm: total normalization
x_mean: expected x value
x2_mean: expected x^2 value
*/
void printLineOnFile(FILE *f, int N, double time, double _Complex *psi, double norm, double x_mean, double x2_mean);

/*
Prints line on file
f: file pointer
N: wave function array length
psi_print: 0 if no, 1 if yes
*/
void printHeaderOnFile(FILE *f, int N, int psi_print);

/*
Executes a C-N step
*/
void crankNicolsonStep(int N, double complex **A_half, double complex *psi, double complex *alpha, double complex *beta, double complex gamma);

/*
Executes a T-S step
*/
void trotterSuzukiStep(int N, double *V, double dt, double complex *psi, double complex eta);

#endif