#include"minimization.h"

double F(gsl_vector * x){
	double s = gsl_vector_get(x,0);
	return s*s;
}

void print_simplex(gsl_matrix * S){
	int n=S->size1,m=n+1;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		printf("%.3g ",gsl_matrix_get(S,i,j));
		}
		printf("\n");
	}	
}

int main(){
	gsl_matrix * simplex = gsl_matrix_alloc(1,2);
	gsl_matrix_set(simplex,0,0,1);
	gsl_matrix_set(simplex,0,1,-2);
	downhill_simplex(F,simplex,0.1);
	print_simplex(simplex);
}
