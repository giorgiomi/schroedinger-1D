#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <complex.h>

/*
Computes the normalization of the wave function at a fixed time
N: array length
psi: wave function array
dx: spatial step size
*/
double normSquared(int N, double _Complex *psi, double dx);

/*
Computes the normalization of the wave function at a fixed time USING SIMPSON
N: array length
psi: wave function array
dx: spatial step size
*/
double normSquaredSimpson(int N, double _Complex *psi, double dx);

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
*/
void printLineOnFile(FILE *f, int N, double time, double _Complex *psi, double norm);

#endif