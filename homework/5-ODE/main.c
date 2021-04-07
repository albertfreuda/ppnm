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

void threebody(double t,gsl_vector * y, gsl_vector * dydt){
	//Eq. of motion for 3-body-problem
	//Give sensible names:
	double x1    = gsl_vector_get(y, 0);
	double x1dot = gsl_vector_get(y, 1);
	double y1    = gsl_vector_get(y, 2);
	double y1dot = gsl_vector_get(y, 3);
	double x2    = gsl_vector_get(y, 4);
	double x2dot = gsl_vector_get(y, 5);
	double y2    = gsl_vector_get(y, 6);
	double y2dot = gsl_vector_get(y, 7);
	double x3    = gsl_vector_get(y, 8);
	double x3dot = gsl_vector_get(y, 9);
	double y3    = gsl_vector_get(y,10);
	double y3dot = gsl_vector_get(y,11);
	//Calculate forces:
	double F12x = (x1-x2)/pow(((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)),1.5);
	double F12y = (y1-y2)/pow(((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)),1.5);
	double F23x = (x2-x3)/pow(((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)),1.5);
	double F23y = (y2-y3)/pow(((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)),1.5);
	double F31x = (x3-x1)/pow(((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)),1.5);
	double F31y = (y3-y1)/pow(((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)),1.5);
	//Set dydt
	gsl_vector_set(dydt, 0,     x1dot);//x1
	gsl_vector_set(dydt, 1,-F12x+F31x);//x1dot
	gsl_vector_set(dydt, 2,     y1dot);//y1
	gsl_vector_set(dydt, 3,-F12y+F31y);//y1dot
	gsl_vector_set(dydt, 4,     x2dot);//x2
	gsl_vector_set(dydt, 5,-F23x+F12x);//x2dot
	gsl_vector_set(dydt, 6,     y2dot);//y2
	gsl_vector_set(dydt, 7,-F23y+F12y);//y2dot
	gsl_vector_set(dydt, 8,     x3dot);//x3
	gsl_vector_set(dydt, 9,-F31x+F23x);//x3dot
	gsl_vector_set(dydt,10,     y3dot);//y3
	gsl_vector_set(dydt,11,-F31y+F23y);//y3dot
}

int main(){
	gsl_vector* yt = gsl_vector_alloc(2);
	gsl_vector* yh = gsl_vector_alloc(2);
	gsl_vector* err= gsl_vector_alloc(2);

	gsl_vector_set(yt,0,1);	
	gsl_vector_set(yt,1,0);	

	double h = .1;
	
	driver(harm_osc,0,yt,4*M_PI,yh,err,h,0.1,0.1);
	
	//Make new index for pyxplot:
	printf("\n\n");
	gsl_vector* ya    = gsl_vector_alloc(3);
	gsl_vector* yb    = gsl_vector_alloc(3);
	gsl_vector* error = gsl_vector_alloc(3);

	gsl_vector_set(ya,0,5e6);
	gsl_vector_set(ya,1,1e6);
	gsl_vector_set(ya,2,0);

	driver(SIR,0,ya,20,yb,error,h,0.1,0.1);
	
	//Make new index for pyxplot:
	printf("\n\n");
	gsl_vector* ya_3b    = gsl_vector_alloc(12);
	gsl_vector* yb_3b    = gsl_vector_alloc(12);
	gsl_vector* error_3b = gsl_vector_alloc(12);

	//Set initial conditions:
	gsl_vector_set(ya_3b, 0, 0.9700436);//x1
	gsl_vector_set(ya_3b, 1, 0.466203685);//x1dot
	gsl_vector_set(ya_3b, 2,-0.24308753);//y1
	gsl_vector_set(ya_3b, 3, 0.43236573);//y1dot
	gsl_vector_set(ya_3b, 4,-0.9700436);//x2
	gsl_vector_set(ya_3b, 5, 0.466203685);//x2dot
	gsl_vector_set(ya_3b, 6, 0.24308753);//y2
	gsl_vector_set(ya_3b, 7, 0.43236573);//y2dot
	gsl_vector_set(ya_3b, 8, 0);//x3
	gsl_vector_set(ya_3b, 9,-2*0.466203685);//x3dot
	gsl_vector_set(ya_3b,10, 0);//y3
	gsl_vector_set(ya_3b,11,-2*0.43236573);//y3dot

	//Solve ODE:
	driver(threebody,0,ya_3b,20,yb_3b,error_3b,h,0.1,0.1);
	
	//Free memory
	gsl_vector_free(yt);
	gsl_vector_free(yh);
	gsl_vector_free(err);
	gsl_vector_free(ya);
	gsl_vector_free(yb);
	gsl_vector_free(error);
	gsl_vector_free(ya_3b);
	gsl_vector_free(yb_3b);
	gsl_vector_free(error_3b);
	return 0;
}
