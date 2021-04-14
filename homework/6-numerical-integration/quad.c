#include<assert.h>
#include<math.h>

//We sample in four points. Not the midpoint!
//We choose 1/6,2/6,4/6,5/6 because they are easily reused when splitting in half.
//Then we use rectangle and trapezium rule to estimate integral and error.

double quad(double f(double),
		double a,
		double b,
		double f2,
		double f4,
		double abstol,
		double reltol,
		int nrec){
	assert(nrec<10000);
	double f1  = f(a+1*(b-a)/6);
	double f5  = f(a+5*(b-a)/6);
	double q   = (f1+f2+f4+f5)/4*(b-a);
	double Q   = (2*f1+f2+f4+2*f5)/6*(b-a);
	double err = fabs(Q-q);
	double tol = abstol + reltol*fabs(Q);
	if(err<tol) return Q;
	else return quad(f,a,(a+b)/2,f1,f2,abstol/sqrt(2),reltol,nrec+1)
		  + quad(f,(a+b)/2,b,f4,f5,abstol/sqrt(2),reltol,nrec+1);
}

double integral(double f(double),
		double a,
		double b,
		double abstol,
		double reltol){
	int nrec = 0;
	double f2 = f(a+2*(b-a)/6);
	double f4 = f(a+4*(b-a)/6);
	return quad(f,a,b,f2,f4,abstol,reltol,nrec);
}

/*First int(f(x),a,b) -> int(f(y),-1,1) by substituting y = -(a+b)/(a-b)+2x/(a-b)*/
/*Then int(f(x),a,b) = int(f(x(y))(b-a)/2,-1,1)*/
/*After that: y = cos(t) so int(f(x(cos(t)))sin(t)(b-a)/2,0,M_PI)*/
double clenshaw_curtis(double f(double),
		double a,
		double b,
		double abstol,
		double reltol){
	double g(double t){ return f((a+b)/2+(a-b)/2*cos(t))*sin(t)*(b-a)/2;}
	return integral(g,0,M_PI,abstol,reltol);
}
