#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void ngradient(double F(gsl_vector * x),
		gsl_vector * x,
		gsl_vector * grad);

int qnewton(double F(gsl_vector * x), //function to minimize
		gsl_vector * x,        //initial guess
		double eps);            //tolerance
