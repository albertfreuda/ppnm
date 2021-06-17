#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<assert.h>
#include<stdio.h>


void matrix_print(gsl_matrix * M);

void cholesky_decomp(gsl_matrix * A,gsl_matrix * L);

void rand_SPD(gsl_matrix * SPD);
