#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"declarations.h"

double pmc_int(double f(int dim, double * x),
		int dim,
		double * a,
		double * b,
		double * err,
		int N);

double f(int dim, double * x){
	double norm_sq = 0;
	for(int i = 0; i<dim; i++) norm_sq += x[i]*x[i];
	if(norm_sq<=1)	return 1;
	else		return 0;
}

double horrible_function(int dim, double * x){
	return 1./(1-cos(x[0])*cos(x[1])*cos(x[2]))/M_PI/M_PI/M_PI;
}

int main(){
	int runs = 5;
	int N = 1e4;
	double a[3] = {0,0,0}, b[3] = {M_PI,M_PI,M_PI},err_p,err_ld;
	double tru_val = 1.3932039296768591842462603255;
	for (int n=0;n<runs;n++){
	double res  = pmc_int(horrible_function,3,a,b,&err_p,N);
	double res_low_disc  = ldmc_int(horrible_function,3,a,b,&err_ld,N);
	printf("%i %g %g\n",N,err_p,err_ld);
	N *=10;
	}
	return 0;
}
