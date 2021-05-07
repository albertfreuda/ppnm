#include<gsl/gsl_vector.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include"rk_declarations.h"
#include"linsolve_declarations.h"

static int ncalls;

void newton(void f(gsl_vector* x,gsl_vector* fx),
		gsl_vector*x,
		double acc);

void constant_function(gsl_vector* x,gsl_vector* fx){
	gsl_vector_memcpy(fx,x);
}

void non_linear(gsl_vector* x,gsl_vector* fx){
	ncalls++;
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,x1*x1+x2*x2-1);
	gsl_vector_set(fx,1,x1+x2);
}

void Rosenbrock_grad(gsl_vector* x,gsl_vector* fx){
	ncalls++;
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);
	gsl_vector_set(fx,0,2*(1-x1)*(-1)-400*(x2-x1*x1)*x1);
	gsl_vector_set(fx,1,200*(x2-x1*x1));
}

void schrodinger(double r,gsl_vector * y,gsl_vector * dydt,double epsilon){
		gsl_vector_set(dydt,0,gsl_vector_get(y,1));
		double y1 = gsl_vector_get(y,0);
		gsl_vector_set(dydt,1,-2*(epsilon+1./r)*y1);
	}

void hydrogen(gsl_vector* x,gsl_vector* M){
	ncalls++;
	//Where do we evaluate the wavefunction?
	double rmax = 8;
	
	//This will be the guess for the energy
	double epsilon = gsl_vector_get(x,0);

	
	gsl_vector * ya  = gsl_vector_alloc(2);
	gsl_vector * yb  = gsl_vector_alloc(2);
	gsl_vector * err = gsl_vector_alloc(2);
	
	//Since the ODE diverges for r=0, choose small r for initial condition instead:
	double rmin = 1e-3;
	//But what are the initial conditions in this place? We calculate using limiting case:
	double Fmin = rmin - rmin*rmin;
	double Fminprime = 1-2*rmin;

	gsl_vector_set(ya,0,Fmin);
	gsl_vector_set(ya,1,Fminprime);

	double h=0.1,acc = 0.01, eps = 0.01;
	driver(schrodinger,rmin,ya,rmax,yb,err,h,acc,eps,epsilon);
	gsl_vector_set(M,0,gsl_vector_get(yb,0));
}

void hydrogen2(gsl_vector* x,gsl_vector* M){
	ncalls++;
	//Where do we evaluate the wavefunction?
	double rmax = 5;
	
	//This will be the guess for the energy
	double epsilon = gsl_vector_get(x,0);

	
	gsl_vector * ya  = gsl_vector_alloc(2);
	gsl_vector * yb  = gsl_vector_alloc(2);
	gsl_vector * err = gsl_vector_alloc(2);
	
	//Since the ODE diverges for r=0, choose small r for initial condition instead:
	double rmin = 1e-3;
	//But what are the initial conditions in this place? We calculate using limiting case:
	double Fmin = rmin - rmin*rmin;
	double Fminprime = 1-2*rmin;

	gsl_vector_set(ya,0,Fmin);
	gsl_vector_set(ya,1,Fminprime);

	double h=0.1,acc = 0.01, eps = 0.01;
	double k = sqrt(-2*epsilon);
	driver(schrodinger,rmin,ya,rmax,yb,err,h,acc,eps,epsilon);
	gsl_vector_set(M,0,gsl_vector_get(yb,0)-rmax*exp(-k*rmax));
}


void harm_osc(double t,gsl_vector* y, gsl_vector* dydt,double param){
	//Harmonic oscillator
	//(y1,y2)' = (y2,-y1) = [0,1;-1,0](y1,y2)
	gsl_vector_set(dydt,0, gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-param*param*gsl_vector_get(y,0));
}

int main(){
	gsl_vector * x = gsl_vector_alloc(2);
	double acc = .01;

	printf("Welcome to Root Finding!\n\nThis is part A:\n\n");

	double x0 = -2;
	double y0 =  8;

	gsl_vector_set(x,0,x0);
	gsl_vector_set(x,1,y0);	

	ncalls = 0;
	newton(Rosenbrock_grad,x,acc);
	printf("Search for Rosenbrock valley extrema started at x=%g and y=%g\n",x0,y0);
	printf("Extremum for the Rosenbrock valley at:\n");
	for(int i=0;i<2;i++){
		printf("%g\n",gsl_vector_get(x,i));
	}
	printf("after calling the function %i times.\n",ncalls);

	printf("\nNow comes part B\n\n");
	gsl_vector * y = gsl_vector_alloc(1);
	gsl_vector_set(y,0,-2);

	ncalls = 0;

	printf("Search for energy minimum started at epsilon = %g\n",gsl_vector_get(y,0));
	newton(hydrogen,y,acc);
	
	printf("Energy of s-wave is: %g\n",gsl_vector_get(y,0));
	printf("after calling the function %i times.\n",ncalls);

	//We also need to plot the resulting wave function, so we solve SE again
	double rmax = 5;
	double epsilon = gsl_vector_get(y,0);
	gsl_vector * ya  = gsl_vector_alloc(2);
	gsl_vector * yb  = gsl_vector_alloc(2);
	gsl_vector * err = gsl_vector_alloc(2);
	double rmin = 1e-3;
	double Fmin = rmin - rmin*rmin;
	double Fminprime = 1-2*rmin;
	gsl_vector_set(ya,0,Fmin);
	gsl_vector_set(ya,1,Fminprime);
	double h=0.1, eps = 0.01;
	driver2(schrodinger,rmin,ya,rmax,yb,err,h,acc,eps,epsilon);
	
	printf("\nNow comes part C\n\n");
	printf("Because we use a better boundary condition,\nwe may use a smaller rmax.\n");

	gsl_vector_set(y,0,-2);
	ncalls = 0;
	printf("Search for energy minimum started at epsilon = %g\n",gsl_vector_get(y,0));
	newton(hydrogen2,y,acc);
	
	printf("Energy of s-wave is: %g\n",gsl_vector_get(y,0));
	printf("after calling the function %i times.\n",ncalls);
	
	gsl_vector_free(x);
	gsl_vector_free(y);
	return 0;
}
