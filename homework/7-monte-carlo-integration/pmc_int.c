#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#define RND ((double)rand()/RAND_MAX)

double pmc_int(double f(int dim, double * x),
		int dim,
		double * a,
		double * b,
		double * err,
		int N){
/* To perform MC integration we sample our function in random places*/
// Make volume right size:
double V = 1; for(int i=0;i<dim;i++) V*=b[i]-a[i];
double x[dim],sum=0,sum2=0;
for(int i=0;i<N;i++){
	// Select point:
	for(int j=0;j<dim;j++) {
		x[j] = a[j] + RND*(b[j]-a[j]);
	}
	double f_val = f(dim,x);
	sum+=f_val; sum2 += f_val*f_val;
}
double sigma = sqrt(sum2/N-sum*sum/N/N);
*err = V*sigma/sqrt(N);
return V*sum/N;

}

