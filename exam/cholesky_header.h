#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<assert.h>
#include<stdio.h>

void vector_print(gsl_vector * v);

void lineq_print(gsl_matrix * M, gsl_vector * v);

void matrix_print(gsl_matrix * M);

void cholesky_decomp(gsl_matrix * A);

void rand_SPD(gsl_matrix * SPD);

double cholesky_det(gsl_matrix * A);

void cholesky_linsolve(gsl_matrix * A, gsl_vector * b, gsl_vector * x);

void cholesky_inverse(gsl_matrix * A, gsl_matrix * B);
