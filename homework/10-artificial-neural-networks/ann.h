#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

typedef struct {int n;double (*f)(double);gsl_vector * params;} ANN;

ANN* ANN_alloc(int n,double (*f)(double));

void ANN_free(ANN* network);

double ANN_response(double x,ANN* network);

void ANN_learn(ANN* network, 
		gsl_vector* xs, 
		gsl_vector* ys);
