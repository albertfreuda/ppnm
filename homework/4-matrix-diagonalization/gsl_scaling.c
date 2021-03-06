#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

void timesJ(gsl_matrix* A, int p, int q, double theta);
void Jtimes(gsl_matrix* A, int p, int q, double theta);

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

void print_vector(gsl_vector* x){
	for(int i=0;i<x->size;i++) printf("%g\n",gsl_vector_get(x,i));
}

void jacobi_diag(gsl_matrix* A,gsl_matrix* V);

int main(int argc,char** argv){
	// Specify default size:
	int N=4,M=N;
	// If argument is supplied, use it:
	if(argc>1) N = (int)atof(argv[1]),M=N;
	//Allocate memory:
	gsl_matrix * A = gsl_matrix_alloc(N,M);
	gsl_matrix * E = gsl_matrix_alloc(N,M);
	gsl_vector * v = gsl_vector_alloc(N);
	//Fill matrix A randomly but symmetrically:
	for(int i=0;i<N;i++){
		for(int j=i;j<M;j++){
		double A_ij = 100*(double)rand()/RAND_MAX;
		gsl_matrix_set(A,i,j,A_ij);
		gsl_matrix_set(A,j,i,A_ij);
		}
	}
	// Diagonalize A:
	gsl_eigen_symmv_workspace * ws = gsl_eigen_symmv_alloc(N);
	gsl_eigen_symmv(A,v,E,ws);

	// Free memory:
	gsl_matrix_free(A);
	gsl_matrix_free(E);
	gsl_vector_free(v);
	gsl_eigen_symmv_free(ws);
	return 0;
}
