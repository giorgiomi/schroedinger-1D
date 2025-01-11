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
Computes the normalization of a function at a fixed time USING SIMPSON
N: array length
f: function array
dx: spatial step size
*/
double normSquaredSimpson(int N, double *f, double dx);

/*
Computes the normalization of a function at a fixed time USING TRAPEZOIDS
N: array length
f: function array
dx: spatial step size
*/
double normSquaredTrap(int N, double *f, double dx);

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

#endif