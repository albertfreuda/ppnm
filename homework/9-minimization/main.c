#include"minimization.h"

static int ncalls;
double F(gsl_vector * x){
	ncalls++;
	double s = gsl_vector_get(x,0);
	double t = gsl_vector_get(x,1);
	return s*s+t*t;
}

double Rosenbrock(gsl_vector * v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}

double Himmelblau(gsl_vector * v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}

void minimum_print(const char * title,
		gsl_vector * guess,
		gsl_vector * result,
		int nsteps,
		double tol,
		double val){
	printf("%s\n",title);
	
	printf("Initial guess: ( ");
	for(int i=0;i<guess->size;i++) printf("%.3g ",gsl_vector_get(guess,i));
	printf(")\n");
	
	printf("Found minimum: ( ");
	for(int i=0;i<result->size;i++) printf("%.3g ",gsl_vector_get(result,i));
	printf(")\n");

	printf("Steps        : %i\n",nsteps);
	printf("Accuracy     : %g\n",tol);
	printf("Value at min : %g\n",val);
}

double BreitWigner(gsl_vector * v,double E){
	double m = gsl_vector_get(v,0);
	double Gamma = gsl_vector_get(v,1);
	double A = gsl_vector_get(v,2);
	return A/((E-m)*(E-m)+Gamma*Gamma/4);
}

int main(){
	ncalls = 0;
	gsl_vector * x  = gsl_vector_alloc(2);
	gsl_vector * xc = gsl_vector_alloc(2);

	double tol = 1e-3;

	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,1);

	gsl_vector_memcpy(xc,x);

	int nsteps = qnewton(F,xc,tol);	
	minimum_print("EUCLIDIAN NORM FUNCTION",x,xc,nsteps,tol,F(xc));

	gsl_vector_set(x,0,3);
	gsl_vector_set(x,1,3);
	gsl_vector_memcpy(xc,x);

	nsteps = qnewton(Rosenbrock,xc,tol);	
	minimum_print("ROSENBROCK'S FUNCTION",x,xc,nsteps,tol,Rosenbrock(xc));
	
	gsl_vector_set(x,0,10);
	gsl_vector_set(x,1,10);
	gsl_vector_memcpy(xc,x);

	nsteps = qnewton(Himmelblau,xc,tol);	
	minimum_print("HIMMELBLAU'S FUNCTION",x,xc,nsteps,tol,Himmelblau(xc));

	printf("\nNow comes part B of the exercise:\n");
	//Data
	double Ecm[30] = {101,103,105,107,109,111,113,115,117,119,121,123,125,127,129,131,133,135,137,139,141,143,145,147,149,151,153,155,157,159};
	double sigma[30] = {-0.25,-0.3,-0.15,-1.71,0.81,0.65,-0.91,0.91,0.96,-2.52,-1.01,2.01,4.83,4.58,1.26,1.01,-1.26,0.45,0.15,-0.91,-0.81,-1.41,1.36,0.50,-0.45,1.61,-2.21,-1.86,1.76,-0.50};
	double err[30] = {2,2,1.9,1.9,1.9,1.9,1.9,1.9,1.6,1.6,1.6,1.6,1.6,1.6,1.3,1.3,1.3,1.3,1.3,1.3,1.1,1.1,1.1,1.1,1.1,1.1,1.1,0.9,0.9,0.9};

	//Parameter vector and initial guess:
	gsl_vector * params  = gsl_vector_alloc(3);
	gsl_vector * paramsc = gsl_vector_alloc(3);
	gsl_vector_set(params,0,126);
	gsl_vector_set(params,1,4);
	gsl_vector_set(params,2,20);

	double deviation(gsl_vector * v){
		double res=0;
		for(int i=0;i<30;i++){
res+=(BreitWigner(v,Ecm[i])-sigma[i])*(BreitWigner(v,Ecm[i])-sigma[i])/err[i]/err[i];
		}
		return res;
	}

	nsteps = qnewton(deviation,paramsc,tol);
	minimum_print("HIGGS BOSON FIT",params,paramsc,nsteps,tol,deviation(paramsc));

	//Make index for pypxlot:
	printf("\n\n");
	for(int i=0;i<30;i++){
	printf("%g %g %g %g\n",Ecm[i],sigma[i],BreitWigner(paramsc,Ecm[i]),BreitWigner(params,Ecm[i]));
	}

	return 0;
}
