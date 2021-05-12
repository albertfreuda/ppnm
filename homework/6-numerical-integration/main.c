#include<stdio.h>
#include"quad_declarations.h"
#include<math.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>


double f(double x){
	return sqrt(x);
}

double g(double x){
	return 4*sqrt(1-x*x);
}

double h(double x){
	return 1./sqrt(x);
}

int main(){
	double abstol = 0.01, reltol = 0.01;
	int n_evals;
	double Q;
	Q = integral(h,0,1,abstol,reltol,&n_evals);
	printf("1st part of exercise B: Comparison between number of evaluations\n");
printf("Integral on (0,1) of 1/sqrt(x)=%g and integrand was evaluated %i times\n",Q,n_evals);
	Q = clenshaw_curtis(h,0,1,abstol,reltol,&n_evals);
printf("Integral on (0,1) of 1/sqrt(x)=%g and integrand was evaluated %i times using the Clenshaw-Curtis transformation.\n",Q,n_evals);
	printf("We see that the Clenshaw-Curtis transformation reduces the number of evaluations.\n");
	printf("To further investigate this, see figure out.evaluations.png\n");
	int calls;
	//Change index for plotting INDEX 1.
	printf("Now comes some data for plotting the figure mentioned above\n");
	printf("\n\n");
	//Nested function (only gcc supports this, outside base c)
	double pi_fun(double x){calls++;return 4*sqrt(1-x*x);}
	
	reltol = 0;//Focus on abstol
	double Q_CC;

	//First we calculate #of evaluations for standard integral function:
	for(int i = 0; i<10; i++){
		calls = 0;
		Q = integral(pi_fun,0,1,abstol,reltol,&n_evals);
		printf("%g %i\n",abstol,calls);
		abstol/=10;
	}
	printf("\n\n");//Change index for pyxplot INDEX 2
	//Then we use Clenshaw Curtis integration routine:
	abstol = 0.01;
	for(int i = 0; i<10; i++){
		calls = 0;
		Q_CC = clenshaw_curtis(pi_fun,0,1,abstol,reltol,&n_evals);
		printf("%g %i\n",abstol,calls);
		abstol/=10;
	}
	printf("\n\n");//Change index for pyxplot INDEX 3
	//Now we use GSL:	
	abstol = 0.01;

int size = 100000;
double pi_fun2(double x,void * params){calls++;return 4*sqrt(1-x*x);}
gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(size);
gsl_function pi_gsl;
pi_gsl.function = &pi_fun2;
double Q_gsl;
double err_gsl;

	for(int i = 0; i<12; i++){
		calls = 0;
gsl_integration_qags(&pi_gsl,0,1,abstol,reltol,size,workspace,&Q_gsl,&err_gsl);
		printf("%g %i\n",abstol,calls);
		abstol/=10;
	}
	printf("\n\n");//index for pyxplot INDEX 4
	printf("Here we print the value of the integral using first the ordinary quadratures, then Clenshaw-Curtis and lastly GSL:\n");
	printf("%.25g %.25g %.25g\n",Q,Q_CC,Q_gsl);
	printf("\n\n");//index for pyxplot INDEX 5

	printf("C: Testing integrals with infinite limits (and the error estimates)\n");

	calls = 0;	
	abstol = 0.01;
	reltol = 0.01;
	double b = INFINITY;
	double b_inf(double x){calls++;return 1.0/(x*x);}
	double err;
	Q = integralC(b_inf,1,b,abstol,reltol,&err);	
	printf("Integral of 1/x^2 from 1 to infinity          = %g +- %g (%i evals)\n",Q,err,calls);

	calls = 0;
	double b_inf2(double x,void* param){calls++;return 1.0/(x*x);}
	pi_gsl.function = &b_inf2;
	gsl_integration_qagiu(&pi_gsl,1,abstol,reltol,size,workspace,&Q_gsl,&err_gsl);
	printf("Integral of 1/x^2 from 1 to infinity (gsl)    = %g +- %g (%i evals)\n",Q_gsl,err_gsl,calls);

	calls = 0;
	double gaussian(double x){calls++;return exp(-x*x);}
	Q = integralC(gaussian,-INFINITY,INFINITY,abstol,reltol,&err);
	printf("Integral of gaussian from -inf to inf         = %g +- %g (%i evals)\n",Q,err,calls);

	calls = 0;
	double gaussian2(double x,void* param){calls++;return exp(-x*x);}
	pi_gsl.function = &gaussian2;
	gsl_integration_qagi(&pi_gsl,abstol,reltol,size,workspace,&Q_gsl,&err_gsl);
	printf("Integral of gaussian from -inf to inf (gsl)   = %g +- %g (%i evals)\n",Q_gsl,err_gsl,calls);
	
	calls = 0;
	double a_inf(double x){calls++;return exp(1.0/x)/x/x;}
	Q = integralC(a_inf,-INFINITY,0,abstol,reltol,&err);
	printf("Integral of exp(1/x)/x^2 from -inf to 0       = %g +- %g (%i evals)\n",Q,err,calls);

	calls = 0;
	double a_inf2(double x,void* param){calls++;return exp(1.0/x)/x/x;}
	pi_gsl.function = &a_inf2;
	gsl_integration_qagil(&pi_gsl,0,abstol,reltol,size,workspace,&Q_gsl,&err_gsl);
	printf("Integral of exp(1/x)/x^2 from -inf to 0 (gsl) = %g +- %g (%i evals)\n",Q_gsl,err_gsl,calls);
return 0;
}
