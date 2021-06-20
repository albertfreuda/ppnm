#include"cholesky_header.h"

int main(){
printf("Hello and welcome to my exam project on Cholesky decomposition.\n");
printf("I have implemented the Cholesky-Banachiewicz algorithm.\n");
printf("It performs the decomposition 'in-place'.\n\n");

	//Allocate memory for two matrices
	gsl_matrix * A = gsl_matrix_alloc(3,3);
	gsl_matrix * L = gsl_matrix_alloc(3,3);
	gsl_matrix * T = gsl_matrix_alloc(3,3);
	
	//Generate symmetric positive definite real matrix
	rand_SPD(A);

	//Print it out, to see
	printf("First I generate a symmetric positive definite matrix:\n");
	matrix_print(A);
	gsl_matrix_memcpy(L,A);

	//Perform Cholesky decomposition and print result
	printf("I then perform the Cholesky decomposition. The result L is:\n");
	cholesky_decomp(L);
	matrix_print(L);

	//Test if decomposition worked
	printf("To test if the decomposition worked, I calculate LL':\n");
	gsl_matrix_memcpy(T,L);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,L,T,0,A);

	matrix_print(A);
	printf("which is clearly equal to the matrix we started with.\n");
	gsl_matrix_memcpy(L,A);

	printf("I also calculate the determinant of the matrix:\n");
	double detA = cholesky_det(L);
	printf("Determinant of A is %g\n\n",detA);

	printf("We might also want to solve a linear equation system.\n");
	gsl_vector * b = gsl_vector_alloc(3);
	gsl_vector * x = gsl_vector_alloc(3);

	for(int i = 0; i < b->size; i++){
		gsl_vector_set(b,i,i+1);
	}

	lineq_print(A,b);
	gsl_matrix_memcpy(L,A);

	cholesky_linsolve(L,b,x);

	printf("The solution is found to be:\n");
	vector_print(x);

	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,b);
	
	printf("while the product Ax (that should be equal to %g,%g,%g) is:\n",gsl_vector_get(b,0),gsl_vector_get(b,1),gsl_vector_get(b,2));
	vector_print(b);

	printf("I also inverted the matrix:\n");
	gsl_matrix_memcpy(T,A);
	cholesky_inverse(T,L);

	matrix_print(L);

	printf("To test whether this worked as intented, I calculate A*A^(-1):\n");

	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,A,L,0,T);

	matrix_print(T);

	printf("This is the identity, indicating that everything worked as intented.\n");

	gsl_matrix_free(A);
	gsl_matrix_free(L);
	gsl_matrix_free(T);
	gsl_vector_free(b);
	gsl_vector_free(x);
}
