#include"ann.h"


double activation(double x){
	return x*exp(-x*x);
}
double activation_int(double x){
	return -0.5*exp(-x*x);
}
double activation_deriv(double x){
	return exp(-x*x)-2*x*x*exp(-x*x);
}

double some_function(double x){
	return sin(5*x-1)*exp(-x*x);
}
double some_function_deriv(double x){
	return -exp(-x*x)*(5*sin(5*x-1)+2*x*cos(5*x-1));
}
double some_function_int(double x){
	return cos(5*x-1)*exp(-x*x);
}

int main(){
	//Initialize the neural network:
	//How many deep neurons?
	int n = 4;
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
	
	//Set initial values for parameters:
	for(int i=0;i<my_network->n;i++){
	//Assign a
	gsl_vector_set(my_network->params,3*i+0,a+i*(b-a)/(my_network->n-1));
	//Assign b
	gsl_vector_set(my_network->params,3*i+1,0.5);
	//Assign w
	gsl_vector_set(my_network->params,3*i+2,0.5);
	}
	/*printf("\n\n");
	printf("a0 = %g\n ",gsl_vector_get(my_network->params,0));
	printf("b0 = %g\n ",gsl_vector_get(my_network->params,1));
	printf("w0 = %g\n ",gsl_vector_get(my_network->params,2));
	printf("a1 = %g\n ",gsl_vector_get(my_network->params,3));
	printf("b1 = %g\n ",gsl_vector_get(my_network->params,4));
	printf("w1 = %g\n ",gsl_vector_get(my_network->params,5));
	printf("a2 = %g\n ",gsl_vector_get(my_network->params,6));
	printf("b2 = %g\n ",gsl_vector_get(my_network->params,7));
	printf("w2 = %g\n ",gsl_vector_get(my_network->params,8));
	printf("a3 = %g\n ",gsl_vector_get(my_network->params,9));
	printf("b3 = %g\n ",gsl_vector_get(my_network->params,10));
	printf("w3 = %g\n ",gsl_vector_get(my_network->params,11));
	*//*gsl_vector_set(my_network->params,0,-1.3);
	gsl_vector_set(my_network->params,1,0.68);
	gsl_vector_set(my_network->params,2,1.8);
	gsl_vector_set(my_network->params,3,-0.09);
	gsl_vector_set(my_network->params,4,.49);
	gsl_vector_set(my_network->params,5,2.9);
	gsl_vector_set(my_network->params,6,-0.34);
	gsl_vector_set(my_network->params,7,0.56);
	gsl_vector_set(my_network->params,8,0.22);
	gsl_vector_set(my_network->params,9,1.07);
	gsl_vector_set(my_network->params,10,0.61);
	gsl_vector_set(my_network->params,11,2.04);
*/
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
		double Fval  = ANN_response(xval,my_network);
		double integ = ANN_eval_int(my_network,activation_int,xval);
		double deriv = ANN_eval_deriv(my_network,activation_deriv,xval);
		double trued = some_function_deriv(xval);
		printf("%g %g %g %g %g\n",xval,Fval,integ,deriv,trued);
		xval+=(b-a)/(N-1);
	}
	/*printf("\n\n");
	printf("a0 = %g\n ",gsl_vector_get(my_network->params,0));
	printf("b0 = %g\n ",gsl_vector_get(my_network->params,1));
	printf("w0 = %g\n ",gsl_vector_get(my_network->params,2));
	printf("a1 = %g\n ",gsl_vector_get(my_network->params,3));
	printf("b1 = %g\n ",gsl_vector_get(my_network->params,4));
	printf("w1 = %g\n ",gsl_vector_get(my_network->params,5));
	printf("a2 = %g\n ",gsl_vector_get(my_network->params,6));
	printf("b2 = %g\n ",gsl_vector_get(my_network->params,7));
	printf("w2 = %g\n ",gsl_vector_get(my_network->params,8));
	printf("a3 = %g\n ",gsl_vector_get(my_network->params,9));
	printf("b3 = %g\n ",gsl_vector_get(my_network->params,10));
	printf("w3 = %g\n ",gsl_vector_get(my_network->params,11));
*/
	//Free the allocated memory:
	ANN_free(my_network);	
}
