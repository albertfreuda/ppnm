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
	gsl_matrix * I = gsl_matrix_calloc(N,M);
	gsl_matrix * Acp = gsl_matrix_calloc(N,M);
	//Fill in the A matrix with random numbers:	
	for(int i=0;i<N;i++){
		for(int k=0;k<M;k++){
			gsl_matrix_set(A,i,k,100*(double)rand()/RAND_MAX);
		}
	}
	gsl_matrix_memcpy(Acp,A);
	printf("This is the matrix A which we are trying to decompose:\n");
	print_matrix(A);

	GS_decomp(A,R);
	printf("\nAfter decomposition we have the matrix Q:\n");
	print_matrix(A);//Note that now A->Q and no longer contain og. A
	printf("\nand the matrix R:\n");
	print_matrix(R);//while R is actually R.
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,I);
	printf("\nR is definately upper triangular. \n\nWe can test our decomposition by computing Q^T*Q which should be the identity:\n");
	print_matrix(I);//A is actually Q now, remember? And I=Q'*Q

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,I);
	printf("\nWe might also want to know if QR = A:\n");
	print_matrix(I);//But now I=QR=A
	printf("\nWe conclude that our algorithm works.\n\n");
	printf("Now we solve a system of linear equations. The right hand side is:\n");
	gsl_vector* b = gsl_vector_alloc(N);//These are empty
	gsl_vector* x = gsl_vector_alloc(N);
	for(int i=0;i<N;i++){
		gsl_vector_set(b,i,100*(double)rand()/RAND_MAX);
	}
	vector_print(b);
	printf("\nTo solve it we use our function GS_solve:\n");
	GS_solve(A,R,b,x);//Here A is Q, R is actually R, b is b and solution in x
	printf("\nThe solution is:\n");
	vector_print(x);
	printf("\nTo test our solution we calculate Ax=b:\n");
	gsl_blas_dgemv(CblasNoTrans,1,Acp,x,0,b);
	vector_print(b);//So now Acp is a copy of A (the real A). b is rewritten.
	printf("\nIt seems that our solution works!\n\n");

	printf("Part B:\nFind the the inverse of A.\n");
	gsl_matrix* B = gsl_matrix_alloc(N,M);
	GS_inverse(A,R,B);//A is Q, so takes QR and saves to B
	print_matrix(B);
	printf("\nTo check that we are right, we compute AB=I\n");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Acp,B,0,I);
	print_matrix(I);
	printf("\nand BA=I\n");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,Acp,0,I);
	print_matrix(I);


	gsl_matrix_free(A);gsl_matrix_free(R);gsl_matrix_free(B);
	gsl_matrix_free(I);gsl_matrix_free(Acp);
	gsl_vector_free(b);gsl_vector_free(x);
return 0;
}
