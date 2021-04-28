#include<assert.h>
#include<math.h>
#include<stdio.h>

double corput(int n, int base){
	// Make room for result
	double q = 0, bk = (double) 1./base;
	while(n>0){
		//Add the inverse of digits in backwards order
		q += (n % base)*bk;
		n /= base;
		bk /= base;
	}
	// Return the number in van der Corput sequence
	return q;
}

void halton(int n, int d, double * x){
	// Find prime numbers for bases:
	int base[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
	// Check that we have enough prime numbers for supplying bases.
	int maxd = sizeof(base)/sizeof(int); assert(d <= maxd);
	// Call corput on each entry in x with new prime base:
	for(int i = 0;i<d;i++){
	x[i] = corput(n,base[i]);
	}
}	

void halton2(int n, int d, double * x){
	// Find prime numbers for bases:
	int base[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
	// Check that we have enough prime numbers for supplying bases.
	int maxd = sizeof(base)/sizeof(int); assert(d <= maxd);
	// Call corput on each entry in x with new prime base:
	for(int i = 0;i<d;i++){
	x[i] = corput(n,base[i]);
	}
}	

double ldmc_int(double f(int dim, double * x),
		int dim,
		double * a,
		double * b,
		double * err,
		int N){
/* To perform MC integration we sample our function in random places*/
// Make volume right size:
double V = 1; for(int i=0;i<dim;i++) V*=b[i]-a[i];
double x[dim],y[dim],sum=0,sum2=0;
int realN = N/2;
for(int i=0;i<realN;i++){
	// Select point using halton-corput sequence. Make room
	double RND[dim];
	double RND2[dim];
	// Call ith numbers in sequences and save to RND array
	halton(i+1,dim,RND);
	halton2(i+1,dim,RND2);
	for(int j=0;j<dim;j++) {
		x[j] = a[j] + RND[j]*(b[j]-a[j]);
		y[j] = a[j] + RND2[j]*(b[j]-a[j]);
	}
	double f_val  = f(dim,x);
	double f_val2 = f(dim,y);
	sum+=f_val; sum2 += f_val2;
}
*err = V*fabs(sum-sum2)/realN;
return V*(sum+sum2)/realN/2;

}
