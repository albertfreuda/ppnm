#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"rk_declarations.h"
#include<math.h>

// The function specifying the differential equation
void harm_osc(double t,gsl_vector* y, gsl_vector* dydt){
	//Harmonic oscillator
	//(y1,y2)' = (y2,-y1) = [0,1;-1,0](y1,y2)
	gsl_vector_set(dydt,0, gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

int main(){
	gsl_vector* yt = gsl_vector_alloc(2);
	gsl_vector* yh = gsl_vector_alloc(2);
	gsl_vector* err= gsl_vector_alloc(2);

	gsl_vector_set(yt,0,1);	
	gsl_vector_set(yt,1,0);	

	double h = .1;
	
	driver(harm_osc,0,yt,0.5*M_PI,yh,err,h,0.1,0.1);

	return 0;
}
