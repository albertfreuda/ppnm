#include<stdio.h>
#include<gsl/gsl_interp.h>

double linterp(int n, double* x, double* y, double z);
double linterp_int(int n, double x[], double y[], double z);

int main(){
	// Make data x²
	int n = 5;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i] = i;y[i] = i*i;
	}
	//Print out data
	for(int i=0;i<n;i++){printf("%g %g\n",x[i],y[i]);}
	printf("\n\n");
	//Make gsl interpolation to compare
	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear,n); //Allocate
	gsl_interp_init(linear,x,y,n); //Build
	
	//Prepare for loop
	int N = 100; //no. of points
	double zmin = 0, zmax = 4, delz = (zmax-zmin)/N; //Linspace
	double z=zmin;
	for(int i=0;i<N;i++){ //Loop over all values
		//Calculate y with own interpolation:
		double y_z = linterp(n,x,y,z);
		//and its integral:
		double integral_z = linterp_int(n,x,y,z);
		//Then do the same with gsl-routine:
		double v = gsl_interp_eval(linear,x,y,z,NULL);
		double integ_v=gsl_interp_eval_integ(linear,x,y,x[0],z,NULL);
		//Print out results along with the proper values.
	printf("%g %g %g %g %g %g %g\n",z,y_z,integral_z,z*z,z*z*z/3,v,integ_v);
		//Increment z:
		z += delz;
	}
	return 0;
}

