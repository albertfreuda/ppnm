#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

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

int main(){
	printf("Welcome to my eigen value decomposition homework.\n");
	int N=3,M=N;//We only deal with symmetric quadratic matrices
	gsl_matrix * A = gsl_matrix_alloc(N,M);
	gsl_matrix * V = gsl_matrix_alloc(N,M);
	gsl_matrix * I = gsl_matrix_alloc(N,M);
	//Fill matrix A
	for(int i=0;i<N;i++){
		for(int j=i;j<M;j++){
		double A_ij = 100*(double)rand()/RAND_MAX;
		gsl_matrix_set(A,i,j,A_ij);
		gsl_matrix_set(A,j,i,A_ij);
		}
	}
	printf("Our matrix A:\n");
	print_matrix(A);
	printf("We multiply it by a Jacobi matrix with angle 0.1:\n");
	timesJ(A,0,1,0.1);
	print_matrix(A);
	printf("We complete the Jacobi transformation:\n");
	Jtimes(A,0,1,0.1);
	print_matrix(A);
	printf("Now we make A diagonal:\n");
	jacobi_diag(A,V);
	print_matrix(A);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,V,V,0,I);
	printf("We check if V'V=I:\n");
	print_matrix(I);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(I);

	return 0;
}
