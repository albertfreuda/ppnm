#include"ann.h"


double activation(double x){
	return x*exp(-x*x);
}

double some_function(double x){
	return cos(5*x-1)*exp(-x*x);
}

int main(){
	//Initialize the neural network:
	//How many deep neurons?
	int n = 5;
	ANN* my_network = ANN_alloc(n,activation);

	//printf("Memory allocation succesful for network\n");

	//Generate data for training:
	int N = 30;
	gsl_vector* xs = gsl_vector_alloc(N);	
	gsl_vector* ys = gsl_vector_alloc(N);	
	double a=-1,b=1,delx=(b-a)/(N-1),x;
	x = a;
	for(int i=0;i<N;i++){
		gsl_vector_set(xs,i,x);
		gsl_vector_set(ys,i,some_function(x));
		x+=delx;
	}
	
	//printf("Data generated succesfully\n");
	
	//Set initial values for parameters all to 1:
	for(int i=0;i<my_network->n;i++){
	gsl_vector_set(my_network->params,3*i+0,a+i*(b-a)/(my_network->n-1));
	gsl_vector_set(my_network->params,3*i+1,1);
	gsl_vector_set(my_network->params,3*i+2,1);
	}
	/*gsl_vector_set(my_network->params,0,-0.8);
	gsl_vector_set(my_network->params,1,1);
	gsl_vector_set(my_network->params,2,1);
	gsl_vector_set(my_network->params,3,-0.8);
	gsl_vector_set(my_network->params,4,1);
	gsl_vector_set(my_network->params,5,1);
	gsl_vector_set(my_network->params,6,-0.8);
	gsl_vector_set(my_network->params,7,1);
	gsl_vector_set(my_network->params,8,1);
	gsl_vector_set(my_network->params,9,-0.8);
	gsl_vector_set(my_network->params,10,1);
	gsl_vector_set(my_network->params,11,1);
	gsl_vector_set(my_network->params,12,1);
	gsl_vector_set(my_network->params,13,1);
	gsl_vector_set(my_network->params,14,1);*/

	//printf("Parameters initiated succesfully\n");
	
	//Then try training the network:
	ANN_learn(my_network,xs,ys);

	//printf("Learning executed succesfully\n");

	//Print out data for plotting
	for(int i=0;i<N;i++){
		double xval = gsl_vector_get(xs,i);
		double yval = gsl_vector_get(ys,i);
		printf("%g %g\n",xval,yval);
	}

	printf("\n\n");

	double xval = a;N=100;
	for(int i=0;i<N;i++){
		double Fval = ANN_response(xval,my_network);
		printf("%g %g\n",xval,Fval);
		xval+=(b-a)/(N-1);
	}


	//Free the allocated memory:
	ANN_free(my_network);	
}
