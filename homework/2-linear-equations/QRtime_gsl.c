#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int main(int argc,char** argv){
	int N = 3,M = N;
	if(argc>1) N=(int)atof(argv[1]),M=N;
	gsl_matrix * A = gsl_matrix_alloc(N,M);
	gsl_vector * tau = gsl_vector_alloc(N);
	//Fill in the A matrix with random numbers:	
	for(int i=0;i<N;i++){
		for(int k=0;k<M;k++){
			gsl_matrix_set(A,i,k,100*(double)rand()/RAND_MAX);
		}
	}
	//Solve system using your own routine
	gsl_linalg_QR_decomp(A,tau);	

	gsl_matrix_free(A);gsl_vector_free(tau);
return 0;
}
