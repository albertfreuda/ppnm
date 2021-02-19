#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>

double Gamma(double);

int main(){
	double xmin = 0.1, xmax = 5;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		printf("%10g %10g %10g %10g\n",x,gamma(x),gsl_sf_lngamma(x),log(Gamma(x)));
		}
return 0;
}
