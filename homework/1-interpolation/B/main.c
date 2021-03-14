#include<stdio.h>
#include<math.h>
#include"functions.h"

int main(){
	// Make data
	int n = 6;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i] = (double)i/(n-1);
		y[i] = exp(-x[i]*x[i]);
	}
	//Now we make the spline:
	qspline *s=qspline_alloc(n,x,y);
	int N = 100;
	double zmin = x[0], zmax = x[n-1], delz = (zmax-zmin)/N;
	double z=zmin;
	for(int i=0;i<N;i++){
		double y_z = qspline_eval(s,z);
		double v_z = qspline_integ(s,z);
		double w_z = qspline_deriv(s,z);
		printf("%g %g %g %g %g %g %g\n",z,y_z,v_z,w_z,exp(-z*z),sqrt(M_PI)*erf(z)/2,-2*z*exp(-z*z));
		z += delz;
	}
	qspline_free(s);
	return 0;
}

