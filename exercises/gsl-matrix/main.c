#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX

void vector_print(char s[], gsl_vector * v){
	printf("%s\n",s); // String to print
	// Print each element in v:
	for(int i=0;i<v->size;i++) printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}

void matrix_print(char s[],gsl_matrix * M,gsl_vector * v){
	printf("%s\n",s);
	for(int i=0;i<M->size1;i++){
		printf("[");
		for(int j=0;j<M->size2;j++){
			double Aij = gsl_matrix_get(M,i,j);
			printf("%5g ",Aij);
		}
		printf("]");
		double vi = gsl_vector_get(v,i);
		printf("x%i [%5g]\n",i,vi);
	}
}

int main(){
	int n=3;
	double RHS[3]={6.23,5.37,2.29};
	double LHS[9]={6.13,-2.9,5.86,8.08,-6.31,-3.89,-4.36,1.00,0.19};
	// Memory allocation
	gsl_matrix * A = gsl_matrix_alloc(n,n);
	gsl_matrix * M = gsl_matrix_alloc(n,n);
	gsl_vector * b = gsl_vector_alloc(n);
	gsl_vector * x = gsl_vector_alloc(n);
	gsl_vector * y = gsl_vector_calloc(n);
	// Fill matrices
	for(int i = 0;i < A->size1; i++)
		for(int j=0; j < A->size2; j++)
		{
		gsl_matrix_set(A,i,j,LHS[i+j]);
		}
	for(int i=0;i< b->size; i++){
		gsl_vector_set(b,i,RHS[i]);
		}
	// Solve system of linear equations
	gsl_matrix_memcpy(M,A); // Copy matrix
	gsl_linalg_HH_solve(M,b,x); // Solve
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y); // Check solution
	// Print out two vectors
	matrix_print("The system is:",A,b);	
	vector_print("Solution x:",x);
	vector_print("Right hand side b: ",b);
	vector_print("A*x should be equal to b:",y);
// Free memory
gsl_matrix_free(A);
gsl_matrix_free(M);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);
return 0;
}
