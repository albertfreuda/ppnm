#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

typedef struct {int n; double x;} nx;

// Define function to be integrated
double f (double x, void * params) {
	nx z = *(nx *) params;
	double g = cos(z.n*x-z.x*sin(x))/M_PI;
	return g;
}

// Now integrate the function: Now it uses the integral representation of gamma
double funA(nx z){
	gsl_function F;
	F.function = &f;
	F.params = (void*)&z;
	int limit=999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double a=0,b=M_PI,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	nx t;
	for(t.x=1.0/8;t.x<=15;t.x+=1.0/8){	
	printf("%10g ",t.x);
		for(t.n=0;t.n<=3;t.n+=1){
		printf("%10g ", funA(t));
	}
	printf("\n");
	}
return 0;
}
