#include <stdlib.h>
#include <assert.h>
#include <math.h>
typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;
/* We want different functions:
 * One for building the spline,
 * one for freeing the memory holding the spline,
 * one for evaluating the spline,
 * one for differentiating the spline
 * and one for integrating the spline 
 * We use the following form of the spline-function:
 * s_i(x) = y_i + b_i(x-x_i) + c_i(x-x_i)² + d_i(x-x_i)³
 */
cspline* cspline_alloc(int n, double *x, double *y){
	//This functions allocates memory for and builds a spline.
	//First we allocate memory corresponding to a spline struct:
	cspline *s = (cspline*)malloc(sizeof(cspline));
	//Then we allocate memory for all the parts of the struct:
	s->b = (double*)malloc((n-1)*sizeof(double));//Coefficients
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));//Data
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;//No. of data points
	//Then we fill in the data from the argument of the function:
	for(int i=0;i<n;i++){
		s->x[i] = x[i];s->y[i] = y[i];
	}
	//Now we calculate the coefficients of the spline:
	//First we allocate memory for the linear terms:
	double p[n-1],h[n-1];
	//Then we calculate them:
	for(int i=0;i<n-1;i++){
		//Delta x
		h[i]=x[i+1]-x[i];
		assert(h[i]>0);
		//Slope
		p[i]=(y[i+1]-y[i])/h[i];
	}
	//The remaining coefficients are found from solving matrix eq.
	//First we build the matrices
	double D[n],Q[n-1],B[n];//Make room
	D[0] = 2;//Define starting point, then fill recursively
	for(int i=0;i<n-2;i++){
		D[i+1]=2*h[i]/h[i+1]+2;
	}
	D[n-1]=2;
	Q[0]=1;//Same for Q
	for(int i=0;i<n-2;i++){
		Q[i+1]=h[i]/h[i+1];
	}//And B:
	for(int i=0;i<n-2;i++){
		B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);
	}//Fix endpoint for B:
	B[0]=3*p[0]; B[n-1]=3*p[n-2];
	//Gauss elimination:
	for(int i=1;i<n;i++){
		D[i]-=Q[i-1]/D[i-1];
		B[i]-=B[i-1]/D[i-1];
	}
	//Lock the last b-coefficient:
	s->b[n-1] = B[n-1]/D[n-1];
	//Then use recursion relation to find remaining b_i's:
	for(int i=n-2;i>=0;i--){
		s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	}
	//Then we can calculate c's and d's:
	for(int i=0;i<n-1;i++){
		s->c[i] = (-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
		s->d[i] = (s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
	}
	//Now, since we have assigned both x's, y's, b's and c's of the spline
	//we return the spline object s:
	return s;
}

void cspline_free(cspline *s){
	//We free all parts of struct, until only n is left. Then we free s.
	free(s->x);free(s->y);free(s->b);free(s->c);free(s->d);free(s);
}

double cspline_eval(cspline *s,double z){
	//Check if z is within supported region
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	//Perform binary search to find proper interval:
	int i=0,j=s->n-1;
	while(j-i>1){int m=(i+j)/2; if(z>s->x[m]) i=m; else j=m;}
	double h=z-s->x[i]; //Change variable to relative system
	// Now calculate value using coefficients from spline
	return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}

double cspline_integ(cspline *s,double z){
	double *x=s->x,*y=s->y,*b=s->b,*c=s->c,*d=s->d;
	int n=s->n;
	//Check if z is in supported interval
	assert (z>=x[0] && z<=x[n-1]);
	//Perform binary search for proper interval:
	int j=0,k=s->n-1;
	while(k-j>1){int m=(j+k)/2; if(z>x[m]) j=m; else k=m;}

	//Do integration
	double integral=0;
	for(int i=0;i<k;i++){
double lhs = x[i]*(y[i]-b[i]*x[i]+c[i]*x[i]*x[i]-d[i]*pow(x[i],3)+x[i]*(b[i]/2-c[i]*x[i]+d[i]*3*pow(x[i],2)/2+x[i]*(c[i]/3-d[i]*x[i]+d[i]*x[i]/4)));
	if(z>=x[i+1]){//If we have not yet reached the last interval
double rhs = x[i+1]*(y[i]-b[i]*x[i]+c[i]*x[i]*x[i]-d[i]*pow(x[i],3)+x[i+1]*(b[i]/2-c[i]*x[i]+d[i]*3*pow(x[i],2)/2+x[i+1]*(c[i]/3-d[i]*x[i]+d[i]*x[i+1]/4)));
	integral += rhs - lhs;
	}
	else{
double rhs = z*(y[i]-b[i]*x[i]+c[i]*x[i]*x[i]-d[i]*pow(x[i],3)+z*(b[i]/2-c[i]*x[i]+d[i]*3*pow(x[i],2)/2+z*(c[i]/3-d[i]*x[i]+d[i]*z/4)));
	integral += rhs - lhs;
	}
}

return integral;
}

double cspline_deriv(cspline *s, double z){
	// Pull out for notation:
	double *x=s->x,*b=s->b,*c=s->c,*d=s->d;
	int n=s->n;
	//Check if z is in supported interval
	assert (z>=x[0] && z<=x[n-1]);
	//Perform binary search for proper interval:
	int j=0,k=s->n-1;
	while(k-j>1){int m=(j+k)/2; if(z>x[m]) j=m; else k=m;}
	//Calculate derivative
	double derivative = b[j] + 2*c[j]*(z-x[j])+3*d[j]*pow((z-x[j]),2);	
	return derivative;
}
