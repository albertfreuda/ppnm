#include<stdio.h>
#include"quad_declarations.h"
#include<math.h>
#include<stdio.h>

double f(double x){
	return sqrt(x);
}

double g(double x){
	return 4*sqrt(1-x*x);
}

double h(double x){
	return 1./sqrt(x);
}

int main(){
	double a = 0,b = 1,abstol = 0.001, reltol = 0.01;
	double Q = integral(f,a,b,abstol,reltol);
	printf("The integral of sqrt(x) from 0 to 1 is: %g\n",Q);
	Q = integral(g,a,b,abstol,reltol);
	printf("The integral of 4*sqrt(1-x^2) from 0 to 1 is: %g\n",Q);
	Q = clenshaw_curtis(h,0,1,abstol,reltol);
	printf("The integral of 1/sqrt(x) from 0 to 1 is: %g\n",Q);
	return 0;
}
