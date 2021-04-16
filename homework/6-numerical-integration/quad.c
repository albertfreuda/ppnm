#include<assert.h>
#include<math.h>
#include<stdio.h>

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
		int nrec,
		int* n_evaluations){
	assert(nrec<10000);
	double f1  = f(a+1*(b-a)/6);
	double f5  = f(a+5*(b-a)/6);
	*n_evaluations += 2;
	double q   = (f1+f2+f4+f5)/4*(b-a);
	double Q   = (2*f1+f2+f4+2*f5)/6*(b-a);
	double err = fabs(Q-q);
	double tol = abstol + reltol*fabs(Q);
	if(err<tol) return Q;
	else{
	return quad(f,a,(a+b)/2,f1,f2,abstol/sqrt(2),reltol,nrec+1,n_evaluations)
	     + quad(f,(a+b)/2,b,f4,f5,abstol/sqrt(2),reltol,nrec+1,n_evaluations);
	}
}

double integral(double f(double),
		double a,
		double b,
		double abstol,
		double reltol,
		int* n_evaluations){
	int nrec = 0;
	double f2 = f(a+2*(b-a)/6);
	double f4 = f(a+4*(b-a)/6);
	*n_evaluations = 2;
	double I  = quad(f,a,b,f2,f4,abstol,reltol,nrec,n_evaluations);
	return I;
}

/*First int(f(x),a,b) -> int(f(y),-1,1) by substituting y = -(a+b)/(a-b)+2x/(a-b)*/
/*Then int(f(x),a,b) = int(f(x(y))(b-a)/2,-1,1)*/
/*After that: y = cos(t) so int(f(x(cos(t)))sin(t)(b-a)/2,0,M_PI)*/
double clenshaw_curtis(double f(double),
		double a,
		double b,
		double abstol,
		double reltol,
		int* n_evaluations){
	double g(double t){ return f((a+b)/2+(a-b)/2*cos(t))*sin(t)*(b-a)/2;}
	return integral(g,0,M_PI,abstol,reltol,n_evaluations);
}


double quadC(double f(double),
		double a,
		double b,
		double f2,
		double f4,
		double abstol,
		double reltol,
		int nrec,
		double* err){
	assert(nrec<10000);
	double f1  = f(a+1*(b-a)/6);
	double f5  = f(a+5*(b-a)/6);
	double q   = (f1+f2+f4+f5)/4*(b-a);
	double Q   = (2*f1+f2+f4+2*f5)/6*(b-a);
	*err = fabs(Q-q);
	double tol = abstol + reltol*fabs(Q);
	if(*err<tol) return Q;
	else{
	return quadC(f,a,(a+b)/2,f1,f2,abstol/sqrt(2),reltol,nrec+1,err)
	     + quadC(f,(a+b)/2,b,f4,f5,abstol/sqrt(2),reltol,nrec+1,err);
	}
}


double integralC(double f(double),
		double a,
		double b,
		double abstol,
		double reltol,
		double* err){
	if(isinf(a) == 0){//then a is finite
		if(isinf(b) == 0){//then be is finite 
//printf("a and b are both finite.\n");
int nrec = 0;
double f2 = f(a+2*(b-a)/6);
double f4 = f(a+4*(b-a)/6);
double I  = quadC(f,a,b,f2,f4,abstol,reltol,nrec,err);
return I;
		}
		else{//a is finite and b is infinite
//printf("a is finite and b is infinity.\n");
double g(double t){return f(a+(1-t)/t)/(t*t);}
int nrec = 0;
double g2 = g(2.0/6);
double g4 = g(4.0/6);
double Ig = quadC(g,0,1,g2,g4,abstol,reltol,nrec,err);
return Ig;

		}
	}
	else{//a is -infinity
		if(isinf(b)==0){//b is finite
//printf("a is -infinity and b is finite.\n");
double h(double s){return f(b+s/(s+1))/(1+s)/(1+s);}
int nrec = 0;
double h2 = h(-1+2.0/6);
double h4 = h(-1+4.0/6);
double Ih = quadC(h,-1,0,h2,h4,abstol,reltol,nrec,err);
return Ih;
		}
		else{//a is -infinity and b is infinity
//printf("a is -infinity and b is infinity.\n");
double k(double r){return f(r/(1-r*r))*(1+r*r)/((1-r*r)*(1-r*r));}
int nrec = 0;
double k2 = k(-1+4.0/6);
double k4 = k(-1+8.0/6);
double Ik = quadC(k,-1,1,k2,k4,abstol,reltol,nrec,err);
return Ik;
	}
}
}
