#include"ann.h"
#include"minimization.h"

//Make a function that allocates memory for a network:
ANN* ANN_alloc(int n,double (*f)(double)){
	ANN* network = malloc(sizeof(ANN));
	network->n=n;
	network->f=f;
	network->params = gsl_vector_alloc(3*n);
	return network;
}

void ANN_free(ANN* network){
	gsl_vector_free(network->params);
	free(network);
}

double ANN_response(double x,ANN* network){
	double y=0;
	for(int i=0;i<network->n;i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		y+=network->f((x-a)/b)*w;
	}
	return y;
}

void ANN_learn(ANN* network, 
		gsl_vector* xs, 
		gsl_vector* ys){
	//First make room for a parameter vector p:
	gsl_vector* p = gsl_vector_alloc(network->params->size);
	//Then copy the initial parameter vector into p:
	gsl_vector_memcpy(p,network->params);
	
	//Define cost function:
	double cost(gsl_vector * p){
		//Make sure that the network holds the current parameters:
		gsl_vector_memcpy(network->params,p);
		//Make room for sum:
		double s = 0;
		//Calculate sum of squared deviations:
		for(int k=0;k<xs->size;k++){
		//For each value in xs
		double Fk = ANN_response(gsl_vector_get(xs,k),network);
		double yk = gsl_vector_get(ys,k);
		s+=(Fk-yk)*(Fk-yk);
		}
		return s;
	}

//	printf("The value of the cost function was: %g\n",cost(p));
	/*int N = network->n;	
	gsl_matrix * simplex = gsl_matrix_alloc(3*N,3*N+1);
	for(int i=0;i<3*N+1;i++){
		for(int k=0;k<N;k++){
		gsl_matrix_set(simplex,3*k+0,i,k);
		gsl_matrix_set(simplex,3*k+1,i,1);
		gsl_matrix_set(simplex,3*k+2,i,0.5);
		}
	}*/
	//Minimize the cost function
	double tolerance = 0.001;
	qnewton(cost,p,tolerance);
	//downhill_simplex(cost,simplex,tolerance,p);

//	printf("The value of the cost function is: %g\n",cost(p));

	//Then copy the minimized parameter vector back into network
	gsl_vector_memcpy(network->params,p);
	//and free the temporary vector p
	gsl_vector_free(p);
}

double ANN_eval_int(ANN* network,double F(double x),double x){
	//This function takes in an ANN and a number x.
	//It also takes F - the integral of the activation function
	double sum = 0;
	for(int i=0;i<network->n;i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		sum+=w*b*F((x-a)/b);	
	}
	return sum;
}

double ANN_eval_deriv(ANN* network,double fp(double x),double x){
	//This function takes in an ANN and a number x.
	//It also takes fp - the derivative of the activation function
	double sum = 0;
	for(int i=0;i<network->n;i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		sum+=w*fp((x-a)/b)/b;	
	}
	return sum;
}
