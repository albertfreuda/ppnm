#include<stdio.h>
#include<gsl/gsl_matrix.h>
#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

int print_half_00(gsl_matrix* m){
	double half = 1./2;
	printf( "half m_{00} = %g\n", gsl_matrix_get(m,0,0)*half );
	return 0;
}

int main(void)
{
	gsl_matrix * m = gsl_matrix_alloc(1,1);
	gsl_matrix_set(m,0,0,66);	
	printf("When printed manually m_{00}:%g\n",gsl_matrix_get(m,0,0)/2);
	int status = print_half_00(m);
	if(status>0)
		printf("status=%i : SOMETHING WENT TERRIBLY WRONG (status>0)\n",status);
	else
		printf("status=%i : everything went just fine (status=0)\n",status);
	gsl_matrix_free(m);
return 0;
}
