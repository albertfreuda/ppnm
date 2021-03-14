#include <stdlib.h>
#include <assert.h>
typedef struct {int n; double *x,*y,*b,*c;} qspline;
/* We want different functions:
 * One for building the spline,
 * one for freeing the memory holding the spline,
 * one for evaluating the spline,
 * one for differentiating the spline
 * and one for integrating the spline 
 * We use the following form of the spline-function:
 * s_i(x) = y_i + b_i(x-x_i) + c_i(x-x_i)² 
 * where b_i = p_i - c_i*dx_i */
qspline* qspline_alloc(int n, double *x, double *y){
	//This functions allocates memory for and builds a spline.
	//First we allocate memory corresponding to a spline struct:
	qspline *s = (qspline*)malloc(sizeof(qspline));
	//Then we allocate memory for all the parts of the struct:
	s->b = (double*)malloc((n-1)*sizeof(double));//Coefficients
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));//Data
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;//No. of data points
	//Then we fill in the data from the argument of the function:
	for(int i=0;i<n;i++){
		s->x[i] = x[i];s->y[i] = y[i];
	}
	//Now we calculate the coefficients of the spline:
	//First we allocate memory for the linear terms:
	int i; double p[n-1],h[n-1];
	//Then we calculate them:
	for(i=0;i<(n-1);i++){
		//Delta x
		h[i]=x[i+1]-x[i];
		//Slope
		p[i]=(y[i+1]-y[i])/h[i];
	}
	//Choose first quadratic coefficient arbitrarily (here to 0):
	s->c[0]=0;
	//Then use recursion relation to find others:
	for(i=0;i<n-2;i++){
		s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	}
	//Now we run same recursion but from other way:
	s->c[n-2]/=2; // Take half of the final coefficient
	for(i=n-3;i>=0;i--){
		s->c[i] = (p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	}
	//Then we assign the b-coefficients
	for(i=0;i<n-1;i++){
		s->b[i]=p[i]-s->c[i]*h[i];
	}
	//Now, since we have assigned both x's, y's, b's and c's of the spline
	//we return the spline object s:
	return s;
}

void qspline_free(qspline *s){
	//We free all parts of struct, until only n is left. Then we free s.
	free(s->x);free(s->y);free(s->b);free(s->c);free(s);
}

double qspline_eval(qspline *s,double z){
	//Check if z is within supported region
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	//Perform binary search to find proper interval:
	int i=0,j=s->n-1;
	while(j-i>1){int m=(i+j)/2; if(z>s->x[m]) i=m; else j=m;}
	double h=z-s->x[i]; //Change variable to relative system
	// Now calculate value using coefficients from spline
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

double qspline_integ(qspline *s,double z){
	double *x=s->x,*y=s->y,*b=s->b,*c=s->c;
	int n=s->n;
	//Check if z is in supported interval
	assert (z>=x[0] && z<=x[n-1]);
	//Perform binary search for proper interval:
	int j=0,k=s->n-1;
	while(k-j>1){int m=(j+k)/2; if(z>x[m]) j=m; else k=m;}

	//Do integration
	double integral=0;
	for(int i=0;i<k;i++){
double lhs = x[i]*(y[i]+x[i]*(c[i]*x[i]-b[i]/2-c[i]*2*x[i]/3));
if(z>=x[i+1]){//If we have not yet reached the last interval
double rhs = x[i+1]*(y[i]+x[i]*(c[i]*x[i]-b[i])+x[i+1]*(b[i]/2-c[i]*(x[i]-x[i+1]/3)));
integral += rhs - lhs;
}
else {double rhs = z*(y[i]+x[i]*(c[i]*x[i]-b[i])+z*(b[i]/2-c[i]*(x[i]-z/3)));
integral += rhs - lhs;
}
}

return integral;
}

double qspline_deriv(qspline *s, double z){
	// Pull out for notation:
	double *x=s->x,*b=s->b,*c=s->c;
	int n=s->n;
	//Check if z is in supported interval
	assert (z>=x[0] && z<=x[n-1]);
	//Perform binary search for proper interval:
	int j=0,k=s->n-1;
	while(k-j>1){int m=(j+k)/2; if(z>x[m]) j=m; else k=m;}
	//Calculate derivative
	double derivative = b[j] + 2*c[j]*(z-x[j]);	
	return derivative;
}
