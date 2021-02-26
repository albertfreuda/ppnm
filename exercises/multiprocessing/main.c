#include<pthread.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

typedef struct {int N; unsigned int seed; int N_in;} params;

void* count(void* arg){
	params* args = (params*)arg;
	int N1 = args->N;
	unsigned int seed1 = args->seed;
	int* N_in = args->N_in;
	for(int i=0;i<N1;i++){
		double x = (double)rand_r(&seed1)/(double)RAND_MAX;
		double y = (double)rand_r(&seed1)/(double)RAND_MAX;
		if(pow(x,2)+pow(y,2) <=1){N_in+=1;}
	}
	return NULL;
}

int main(){
	int N = 1e2;
	unsigned int seed1,seed2;
	int N1 = 0,N2 = 0;
	params arg1; arg1.N = N; arg1.seed = seed1; arg1.N_in=N1;
	params arg2; arg1.N = N; arg1.seed = seed2; arg1.N_in=N2;
	//pthread_t t1,t2;
	//pthread_create(&t1,NULL,count,(void*)&arg1);
	//pthread_create(&t2,NULL,count,(void*)&arg2);
	count((void*)arg1);
	double pi = 4*N1/N;
	printf("%g\n",pi);

return 0;
}
