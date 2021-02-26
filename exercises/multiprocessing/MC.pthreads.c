#include<pthread.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(){
	unsigned int seed1;
	int N = 1e6;
	int N_in=0;
	for(int i=0;i<N;i++){
	double x = (double)rand_r(&seed1)/(double)RAND_MAX;
	double y = (double)rand_r(&seed1)/(double)RAND_MAX;
	if(pow(x,2)+pow(y,2) <= 1){N_in+=1;}
	}
	double pi = 4*(double)N_in/(double)N;
	printf("Using N = %i and counting N_in = %i we get pi=%g\n",N,N_in,pi);
return 0;
}
