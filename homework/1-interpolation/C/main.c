#include<stdio.h>
#include<math.h>
#include<gsl/gsl_interp.h>
#include"functions.h"

int main(void){
	// Make data
	int n = 6;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i] = (double)i/(n-1);
		y[i] = exp(-x[i]*x[i]);
	}
	//Now we make the spline:
	cspline *s=cspline_alloc(n,x,y);
	int N = 100;
	//We also make splines from gsl:
	gsl_interp* cspline = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp_init(cspline,x,y,n);
	double zmin = x[0], zmax = x[n-1], delz = (zmax-zmin)/N;
	double z=zmin;
	for(int i=0;i<N;i++){
		//Evaluate our own splines:
		double y_z = cspline_eval(s,z);
		double v_z = cspline_integ(s,z);
		double w_z = cspline_deriv(s,z);
		printf("%g %g %g %g %g %g %g\n",z,y_z,v_z,w_z,exp(-z*z),sqrt(M_PI)*erf(z)/2,-2*z*exp(-z*z));
		z += delz;
	}
	printf("\n\n");
	z=zmin;
	for(int i=0;i<N-1;i++){
		z+=delz;
		double y_z  = gsl_interp_eval(cspline,x,y,z,NULL);
		double integ= gsl_interp_eval_integ(cspline,x,y,x[0],z,NULL);
		double deriv= gsl_interp_eval_deriv(cspline,x,y,z,NULL);
		printf("%g %g %g %g\n",z,y_z,integ,deriv);
	}
	printf("\n\n");
	for(int i=0;i<n;i++){
		printf("%g %g\n",x[i],y[i]);
	}
	cspline_free(s);
	return 0;
}

