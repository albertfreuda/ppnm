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
	//Specify resolution
	int N=40,M=N;//We only deal with symmetric quadratic matrices
	double s = 1.0/(N+1);
	//Allocate memory
	gsl_matrix * H = gsl_matrix_alloc(N,M);
	gsl_matrix * Acp = gsl_matrix_alloc(N,M);
	gsl_matrix * V = gsl_matrix_alloc(N,M);
	gsl_matrix * I = gsl_matrix_alloc(N,M);
	//Fill matrix A
	for(int i=0;i<N-1;i++){
	gsl_matrix_set(H,i,i,-2);
	gsl_matrix_set(H,i,i+1,1);
	gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,N-1,N-1,-2);
	gsl_matrix_scale(H,-1/s/s);
	//Diagonalize matrix:
	jacobi_diag(H,V);
	//Print our eigenvalues:
	for(int k=0;k<15;k++){
		//Calculate the analytic eigenvalue
		double analytic = M_PI*M_PI*(k+1)*(k+1);
		//Find numerical:
		double numeric = gsl_matrix_get(H,k,k);
		printf("%i %g %g \n",k,analytic,numeric);
	}
	//Also, print out eigenfunctions and compare with the expectation:
	printf("\n\n");//To have new index for pyxplot
	printf("0 0 0 0\n");//First point
	for(int i=0;i<N;i++){
		printf("%g %g %g %g\n",(i+1.0)/(N+1),gsl_matrix_get(V,i,0),gsl_matrix_get(V,i,1),gsl_matrix_get(V,i,2));
	}
	printf("1 0 0 0\n");//Last point

	gsl_matrix_free(H);
	gsl_matrix_free(Acp);
	gsl_matrix_free(V);
	gsl_matrix_free(I);

	return 0;
}
