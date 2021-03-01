#include<pthread.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

struct params {unsigned int seed; int N,N_in;};

void* count(void* arg){
	struct params * p  = (struct params*)arg; // Cast argument to proper type
	p->N_in=0; // Set the starting value of N_in to zero
	for(int i=0;i<(p->N);i++){
		double x = (double)rand_r(&(p->seed))/RAND_MAX;
		double y = (double)rand_r(&(p->seed))/RAND_MAX;
		if(x*x + y*y <1)p->N_in++;
	}
	return NULL;
}

int main(int argc, char** argv){ // Optional number of points as argument
	int N = (int)1e2; // Default number of points
	if(argc>1) N=(int)atof(argv[1]); // Define number of points from opt. arg.
	struct params p1 = {.seed = 1,.N=N/3,.N_in=0};
	struct params p2 = {.seed = 2,.N=N/3,.N_in=0};
	struct params p3 = {.seed = 3,.N=N/3,.N_in=0};
	pthread_t t1,t2;
	pthread_create(&t1,NULL,count,(void*)&p1);
	pthread_create(&t2,NULL,count,(void*)&p2);
	count((void*) &p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	int N1 = p3.N_in + p2.N_in + p1.N_in;
	int Ntrue = p1.N + p2.N + p3.N;
	double pi = 4*(double)N1/Ntrue;
	printf("%i %g %g\n",Ntrue,pi,fabs(pi-M_PI));

return 0;
}
