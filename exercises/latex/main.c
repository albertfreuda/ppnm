#include <stdio.h>
#include <math.h>

double ex(double x);

int main(){
	int N = 100;
	double xmin=-5,xmax=5;
	double x=xmin;
	double delx = (xmax - xmin)/N;
	for(int i=0;i<N;i++){
		printf("%g %g %g\n",x,exp(x),ex(x));
		x+=delx;
	}
	return 0;
}
