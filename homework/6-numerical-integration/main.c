#include<stdio.h>
#include"quad_declarations.h"
#include<math.h>
#include<stdio.h>
#include<time.h>


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
	printf("Part of exercise: Comparison between number of evaluations\n");
printf("Integral on (0,1) of 1/sqrt(x)=%g and integrand was evaluated %i times\n",Q,n_evals);
	Q = clenshaw_curtis(h,0,1,abstol,reltol,&n_evals);
printf("Integral on (0,1) of 1/sqrt(x)=%g and integrand was evaluated %i times\n",Q,n_evals);

	//Change index for plotting
	printf("\n\nThe PI approximations.\n");
printf("Tolerance : Value of integral : Time : # of evaluations : Value CC : Time CC : # CC\n");	
	int n_evals_CC;
	abstol =1e-10; 	
	clock_t begin = clock();
	Q = integral(g,0,1,abstol,0,&n_evals);
	clock_t end = clock();
	double time = (double)(end-begin)/CLOCKS_PER_SEC;
	
	begin = clock();
	double Q_CC = clenshaw_curtis(g,0,1,abstol,0,&n_evals_CC);
	end = clock();
	double time_CC = (double)(end-begin)/CLOCKS_PER_SEC;
	printf("Evaluating PI integral without ClenshawCurtis:\n");
	printf("               Q = %.20g\n",Q);	
	printf("            Time = %g seconds\n",time);	
	printf("# of evaluations = %i\n\n",n_evals);	
	printf("Evaluating PI integral with ClenshawCurtis:\n");
	printf("               Q = %.20g\n",Q_CC);	
	printf("            Time = %g seconds\n",time_CC);	
	printf("# of evaluations = %i\n\n",n_evals_CC);	

return 0;
}
