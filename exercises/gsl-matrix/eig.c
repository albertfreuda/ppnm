#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX

void vector_print(char s[], gsl_vector * v){
	printf("%s\n",s); // String to print
	// Print each element in v:
	for(int i=0;i<v->size;i++) printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(char s[],gsl_matrix * M){
	printf("%s\n",s);
	for(int i=0;i<M->size1;i++){
		printf("[");
		for(int j=0;j<M->size2;j++){
			double Aij = gsl_matrix_get(M,i,j);
			printf("%10g ",Aij);
		}
		printf("]");
		printf("\n");
	}
}

int main(){
	int n=3;
	// Memory allocation
	gsl_matrix * A = gsl_matrix_alloc(n,n);
	gsl_matrix * M = gsl_matrix_alloc(n,n);
	gsl_matrix * evec = gsl_matrix_alloc(n,n);
	gsl_vector * x = gsl_vector_alloc(n);
	// Fill matrix
	for(int i = 0;i < A->size1; i++)
		for(int j=0; j < A->size2; j++)
		{
		double Aij = 1.0/(i+j+1);
		gsl_matrix_set(A,i,j,Aij);
		}	
	gsl_matrix_memcpy(M,A); // Copy matrix
	// Solve eigenvalue problem
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
	gsl_eigen_symmv(M,x,evec,w);
	// Print out solution	
	matrix_print("The matrix to diagonalize is:",A);	
	vector_print("Eigenvalues are:",x);
	matrix_print("Eigenvectors are:",evec);
// Free memory
gsl_matrix_free(A);
gsl_matrix_free(evec);
gsl_matrix_free(M);
gsl_vector_free(x);
gsl_eigen_symmv_free(w);
return 0;
}
