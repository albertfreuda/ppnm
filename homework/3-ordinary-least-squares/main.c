#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

double linfun(int i, double x){
	if((i=0)) return 1;
	if((i=1)) return x;
}

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

//void ls_fit(gsl_vector* x,gsl_vector* y, double fitfun(int,double),gsl_vector* sigma,gsl_vector* c){

int main(){
	double xs[9] = {1,2,3,4,6,9,10,13,15};
	double ys[9] = {117,100,88,72,53,29.5,25.2,15.2,11.1};

	gsl_vector * x = gsl_vector_alloc(9);
	gsl_vector * y = gsl_vector_alloc(9);
	gsl_vector * sigma = gsl_vector_alloc(9);
	gsl_vector * c = gsl_vector_alloc(2);
	for(int i=0;i<x->size;i++){
		//printf("%i\n",i);
		gsl_vector_set(x,i,xs[i]);
		gsl_vector_set(y,i,ys[i]);
		gsl_vector_set(sigma,i,ys[i]/20);
	}
	//vector_print(x);
	//vector_print(y);
	//vector_print(sigma);

	//Then we perform the fit
	//ls_fit(x,y,linfun,sigma,c);

	vector_print(c);

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(sigma);
	gsl_vector_free(c);
return 0;
}
