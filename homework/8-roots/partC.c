#include<gsl/gsl_vector.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include"rk_declarations.h"
#include"linsolve_declarations.h"

static int ncalls,rmax;

void newton(void f(gsl_vector* x,gsl_vector* fx),
		gsl_vector*x,
		double acc);

void schrodinger(double r,gsl_vector * y,gsl_vector * dydt,double epsilon){
		gsl_vector_set(dydt,0,gsl_vector_get(y,1));
		double y1 = gsl_vector_get(y,0);
		gsl_vector_set(dydt,1,-2*(epsilon+1./r)*y1);
	}

void hydrogen(gsl_vector* x,gsl_vector* M){
	ncalls++;
	//Where do we evaluate the wavefunction?
	
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


int main(){
	//Make room for the energy
	double acc = 0.01;
	gsl_vector * x  = gsl_vector_alloc(1);
	gsl_vector * x2 = gsl_vector_alloc(1);
	//We want to solve for the energy for several different values of rmax
	for(rmax=2;rmax<10;rmax++){
		gsl_vector_set(x ,0,-3);
		gsl_vector_set(x2,0,-3);
		ncalls = 0;
		newton(hydrogen ,x ,acc);
		//printf("%i %i ",rmax,ncalls);
		ncalls = 0;
		newton(hydrogen2,x2,acc);
		//printf("%i\n",ncalls);
		printf("%i %g %g\n",rmax,gsl_vector_get(x,0),gsl_vector_get(x2,0));
	}
	gsl_vector_free( x);
	gsl_vector_free(x2);
	return 0;
}
