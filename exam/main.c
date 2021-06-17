#include"cholesky_header.h"

int main(){

	gsl_matrix * A = gsl_matrix_alloc(3,3);
	gsl_matrix * L = gsl_matrix_alloc(3,3);
		
	rand_SPD(A);

	matrix_print(A);

	cholesky_decomp(A,L);
	matrix_print(L);

	gsl_matrix_free(A);
	gsl_matrix_free(L);
}
