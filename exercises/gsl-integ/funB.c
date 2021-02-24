#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

// Define function to be integrated
double f (double x, void * params) {
	// For generality we include parameter, though it is not used here
	double z = *(double *) params;
	double g = z*2*exp(-pow(x,2))/sqrt(M_PI);
	return g;
}

// Now integrate the function: Now it uses the integral representation of gamma
double funA(double z){
	gsl_function F;
	F.function = &f;
	double par = 1;
	F.params = (void*)&par;
	int limit=999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double a=0,b=z,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double x=-5.0;x<=5;x+=1.0/8){
		printf("%10g %10g\n",x,funA(x));
	}
return 0;
}
