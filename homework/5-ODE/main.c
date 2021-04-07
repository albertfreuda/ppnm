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

void SIR(double t,gsl_vector* y, gsl_vector* dydt){
	// SIR-model of epidemic development. Time in days
	// y0' = -y1*y0/N*Tc
	// y1' =  y1*y0/N*Tc-y1/Tr
	// y2' =  y1/Tr
	double N  = 6e6; //Population size
	double Tc = 3;
	double Tr = 10;
	
double y0p=-gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc);
double y1p= gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc)-gsl_vector_get(y,1)/Tr;
double y2p= gsl_vector_get(y,1)/Tr;
	
	gsl_vector_set(dydt,0, y0p);
	gsl_vector_set(dydt,1, y1p);
	gsl_vector_set(dydt,2, y2p);
}

int main(){
	gsl_vector* yt = gsl_vector_alloc(2);
	gsl_vector* yh = gsl_vector_alloc(2);
	gsl_vector* err= gsl_vector_alloc(2);

	gsl_vector_set(yt,0,1);	
	gsl_vector_set(yt,1,0);	

	double h = .1;
	
	driver(harm_osc,0,yt,2*M_PI,yh,err,h,0.1,0.1);
	
	//Make new index for pyxplot:
	printf("\n\n");
	gsl_vector* ya    = gsl_vector_alloc(3);
	gsl_vector* yb    = gsl_vector_alloc(3);
	gsl_vector* error = gsl_vector_alloc(3);

	gsl_vector_set(ya,0,5e6);
	gsl_vector_set(ya,1,1e6);
	gsl_vector_set(ya,2,0);

	driver(SIR,0,ya,20,yb,error,h,0.1,0.1);
	
	gsl_vector_free(yt);
	gsl_vector_free(yh);
	gsl_vector_free(err);
	gsl_vector_free(ya);
	gsl_vector_free(yb);
	gsl_vector_free(error);
	return 0;
}
