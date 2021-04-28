#include<gsl/gsl_vector.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include"rk_declarations.h"


void newton(void f(gsl_vector* x,gsl_vector* fx),
		gsl_vector*x,
		double acc);

void constant_function(gsl_vector* x,gsl_vector* fx){
	gsl_vector_memcpy(fx,x);
}

void non_linear(gsl_vector* x,gsl_vector* fx){
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,x1*x1+x2*x2-1);
	gsl_vector_set(fx,1,x1+x2);
}

void Rosenbrock_grad(gsl_vector* x,gsl_vector* fx){
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,2*(1-x1)*(1-x1)-400*(x2-x1*x1)*x1);
	gsl_vector_set(fx,1,200*(x2-x1*x1));
}


void hydrogen(gsl_vector* x,gsl_vector* M){
	double rmax = 8;
	
	double epsilon = gsl_vector_get(x,0);

	void schrodinger(double r,gsl_vector * y,gsl_vector * dydt){
		gsl_vector_set(dydt,0,gsl_vector_get(y,1));
		double y1 = gsl_vector_get(y,0);
		gsl_vector_set(dydt,1,-2*(epsilon+1./r)*y1);
	}
	gsl_vector * ya  = gsl_vector_alloc(2);
	gsl_vector * yb  = gsl_vector_alloc(2);
	gsl_vector * err = gsl_vector_alloc(2);

	gsl_vector_set(ya,0,0);
	gsl_vector_set(ya,1,0);

	double h=0.1,acc = 0.01, eps = 0.01;
	printf("You now start integrating the Schrodinger eq.");
	driver(schrodinger,0,ya,rmax,yb,err,h,acc,eps);
	gsl_vector_set(M,0,gsl_vector_get(yb,0));
	printf("The value at rmax is: %g\n",gsl_vector_get(M,0));
}

int main(){
	gsl_vector * x = gsl_vector_alloc(2);
	double acc = .01;

	gsl_vector_set(x,0,1.5);
	gsl_vector_set(x,1,1.5);	

	newton(Rosenbrock_grad,x,acc);
	printf("Extremum for the Rosenbrock valley at:\n");
	for(int i=0;i<2;i++){
		printf("%g\n",gsl_vector_get(x,i));
	}

	gsl_vector * y = gsl_vector_alloc(1);

	newton(hydrogen,y,acc);
printf("The minimum energy of s-wave was found to be %g\n",gsl_vector_get(y,0));

	gsl_vector_free(x);
	gsl_vector_free(y);
	return 0;
}
