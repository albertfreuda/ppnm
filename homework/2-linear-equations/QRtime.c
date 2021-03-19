#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* A, gsl_matrix* R,gsl_vector* b,gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

void print_matrix(gsl_matrix* A){
	int n=A->size1,m=A->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		double A_ij = gsl_matrix_get(A,i,j);
		if(fabs(A_ij)<1e-12) printf("%12i ",0);
		else printf("%12g ",A_ij);
		}
		printf("\n");
	}
}

void vector_print(gsl_vector* x){
	for(int i=0;i<x->size;i++) printf("%g\n",gsl_vector_get(x,i));
}

int main(int argc,char** argv){
	int N = 3,M = N;
	if(argc>1) N=(int)atof(argv[1]),M=N;
	gsl_matrix * A = gsl_matrix_alloc(N,M);
	gsl_matrix * R = gsl_matrix_alloc(N,M);
	//Fill in the A matrix with random numbers:	
	for(int i=0;i<N;i++){
		for(int k=0;k<M;k++){
			gsl_matrix_set(A,i,k,100*(double)rand()/RAND_MAX);
		}
	}
	GS_decomp(A,R);

	gsl_matrix_free(A);gsl_matrix_free(R);
return 0;
}
