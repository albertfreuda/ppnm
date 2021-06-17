#include"cholesky_header.h"

void matrix_print(gsl_matrix * M){
	int n=M->size1;
	int m=M->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		printf("%g ",gsl_matrix_get(M,i,j));
		}
		printf("\n");
	}
	printf("\n");
}

//This algorithm decomposes the real symmetric positive-definite matrix A into LL^T
void cholesky_decomp(gsl_matrix * A,gsl_matrix * L){
	//Assert that A is at least a square matrix? Just to avoid error...
	
	//Determine dimension of problem
	int n=A->size1;
	//Sum over rows
	for(int i=0;i<n;i++){
		//Sum over entries in row left of diagonal	
		for(int j=0;j<i+1;j++){
			//Calculate the sum in the formula (41) in book
			double sum=0;
			for(int k=0;k<j;k++){
				double Lik = gsl_matrix_get(L,i,k);
				double Lij = gsl_matrix_get(L,i,j);
				sum += Lik*Lij;
			}
			//Special treatment for diagonal elements
			if(i==j){
				double Aii = gsl_matrix_get(A,i,i);
			       	gsl_matrix_set(L,i,j,sqrt(Aii-sum));
			}
			else{
				double Aij = gsl_matrix_get(A,i,j);
				double Ljj = gsl_matrix_get(L,j,j);
				double Lij = 1.0/Ljj*(Aij-sum);
				gsl_matrix_set(L,i,j,Lij);
			}
		}
	}
}

void rand_SPD(gsl_matrix * SPD){
	//This function creates a symmetric positive definite matrix 
	//by generating a diagonal matrix D with only positive entries
	//and a random matrix Q. Then the product Q'DQ is symmetric
	//positive definite.
	int n=SPD->size1;
	int m=SPD->size2;
	//Assert that matrix is square
	assert(n==m);
	
	//Allocate memory for D and Q
	gsl_matrix * D = gsl_matrix_alloc(n,m);
	gsl_matrix * Q = gsl_matrix_alloc(n,m);
	
	//Populate matrices randomly
	for(int i=0;i<n;i++){
		//Set eigenvalues to random positive number
		gsl_matrix_set(D,i,i,10*(double)rand()/RAND_MAX);
		for(int j=0;j<m;j++){
			//Populate random matrix Q
			gsl_matrix_set(Q,i,j,10*(double)rand()/RAND_MAX);
		}
	}
	//Calculate SPD = Q'DQ: First SPD=1*DQ
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,D,Q,0,SPD);
	//Then D = Q'SPD:
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Q,SPD,0,D);
	//Then copy D to SPD:
	gsl_matrix_memcpy(SPD,D);

	//And free temporary matrices:
	gsl_matrix_free(D);
	gsl_matrix_free(Q);
}
