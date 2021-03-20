#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void GS_decomp(gsl_matrix* A,gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

void print_matrix(gsl_matrix* A);

void vector_print(gsl_vector* x);

void ls_fit(gsl_vector* x,gsl_vector* y, double f(int,double),gsl_vector* sigma,gsl_vector* c,gsl_matrix* S){
	//First we make the A-matrix
	int n = x->size, m = c->size;
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_vector* b = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			double Aij = f(j,gsl_vector_get(x,i))/gsl_vector_get(sigma,i);
			gsl_matrix_set(A,i,j,Aij);
		}
		double bi = gsl_vector_get(y,i)/gsl_vector_get(sigma,i);
		gsl_vector_set(b,i,bi);
	}
	//Then we decompose A:
	gsl_matrix* R=gsl_matrix_alloc(m,m);
	gsl_matrix* I=gsl_matrix_alloc(m,m);
	GS_decomp(A,R);
	//and solve the system Rc=Q'b
	GS_solve(A,R,b,c);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,R,R,0,I);
	GS_decomp(I,R);
	GS_inverse(I,R,S);
	//Then the vector c holds all the coefficients.
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(I);
	gsl_vector_free(b);
}
