#include<assert.h>

double linterp_int(int n, double x[], double y[], double z){
	// Check if z is in supported interval
	assert (n>1 && z>=x[0] && z<=x[n-1]);
	int j=0,k=n-1;
	while (k-j>1){int m=(j+k)/2; if(z>x[m]) j=m; else k=m;}	

	double integral=0;
	for(int i=0; i<k; i++){
		// Calculate the slope
		double p = (y[i+1] - y[i]) / (x[i+1] - x[i]);
		// Calculate left edge of integral from linear interpolation.
		double lhs = y[i]*x[i] + p*(x[i]*x[i]/2-x[i]*x[i]);
		if(z>=x[i+1]){
		// Calculate right edge from linear interpolation.
		double rhs = y[i]*x[i+1] + p*(x[i+1]*x[i+1]/2 -x[i]*x[i+1]);
		integral += rhs - lhs;
		}
		else{
		double rhs = y[i]*z + p*(z*z/2 -x[i]*z);
		integral += rhs - lhs;
		}
	}
	return integral;
}
