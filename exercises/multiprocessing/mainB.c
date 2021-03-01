#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<stdlib.h>

int main(int argc, char** args){
	unsigned int x=1,y=2,z=3; // seeds for rand_r
	int N = (int) 4e7; // number of points
	if(argc>1) N=(int)atof(args[1])/3;
	int N1=0,N2=0,N3=0; // memory for hits
#pragma omp parallel sections
	{
	#pragma omp section
		{
		for(int i=0;i<N;i++){
		double a = (double) rand_r(&x)/RAND_MAX;
		double b = (double) rand_r(&x)/RAND_MAX;
		if(a*a+b*b <=1)N1++;
		}
		}
	#pragma omp section
		{
		for(int j=0;j<N;j++){
		double c = (double) rand_r(&y)/RAND_MAX;
		double d = (double) rand_r(&y)/RAND_MAX;
		if(c*c+d*d <=1)N2++;
		}
		}
	#pragma omp section
		{
		for(int k=0;k<N;k++){
		double e = (double) rand_r(&z)/RAND_MAX;
		double f = (double) rand_r(&z)/RAND_MAX;
		if(e*e+f*f <=1)N3++;
		}
		}
	}
	int Nhit = N1 + N2 + N3;
	int Ntrue = 3*N;
	double pi = 4*(double)Nhit/Ntrue;
	printf("%g\n",pi);
return 0;
}
