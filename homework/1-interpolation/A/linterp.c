#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<assert.h>

double linterp(int n, double* x, double* y, double z){
	// Check if z is inside supported region
	assert (n>1 && z>=x[0] && z<=x[n-1]);
	// Initialize indices
	int i=0,j=n-1;
	// Choose indices: Start from the middle m=(n-1)/2.
	// If z is larger than x[m] you set i=m.
	// If z is smaller than x[m] you set j=m.
	// This way, you are moving the index i towards the value were
	// x[i] < z z x[i+1].
	while(j-i>1){int m=(i+j)/2;if(z>x[m]) i=m; else j=m;}
	assert(x[i+1]>x[i]);
	// Calculate slope in interval.
	double p = (y[i+1] - y[i]) / (x[i+1] - x[i]);
	// Calculate z from linear interpolation.
	return y[i] + p * ( z - x[i] );
}
