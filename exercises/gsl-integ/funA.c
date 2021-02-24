#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

// Define your function to be integrated
double f (double x, void * params) {
	// For generality we include parameter, though it is not used here
	double z = *(double *) params;
	double g = z*log(x)/sqrt(x);
	return g;
}

// Now integrate the function: Now it uses the integral representation of gamma
double funA(double z){
	gsl_function F;
	F.function = &f;
	F.params = (void*)&z;
	int limit=999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double a=0,b=1,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	double z = 1;
	printf("The integral of log(x)/sqrt(x) from 0 to 1 is: %g\n",funA(z));	
return 0;
}
