#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector_double.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<stdio.h>
#include<float.h>
#include"linsolve_declarations.h"
#define eps DBL_EPSILON

void newton(void f(gsl_vector* x,gsl_vector* fx),
		gsl_vector*x,
		double acc){
	// Overall structure is the following:
	// First compute jacobian
	// Then solve J delx = -f(x) to find best delx
	// Then optimize the length of the step through backtracking linesearch
	// Then take the step and update your point. 
	// If value of function is larger than tolerance 
	// Recalculate J and start over.
	//
	int dim = x->size;



	//Make room for evaluating function at x, and all delx's
	gsl_vector * xp = gsl_vector_alloc(dim);
	gsl_vector * fx = gsl_vector_alloc(dim);
	gsl_vector * fdx = gsl_vector_alloc(dim);
	gsl_vector * xpdx = gsl_vector_alloc(dim);
	//Make room for Jacobian
	gsl_matrix * J = gsl_matrix_alloc(dim,dim);
	gsl_matrix * Js = gsl_matrix_alloc(dim,dim);
	//Calculate f(x) and store in fx:
	f(x,fx);
	double normfx = gsl_blas_dnrm2(fx);
	while(normfx>acc){
	//Calculate Jacobian matrix:
	//Define step size:
	double delx = sqrt(eps);
	//Now run through matrix:
	for(int k=0;k<dim;k++){
		//Calculate first x+delx in direction k:
		//First we copy the content of x to xp
		gsl_vector_memcpy(xp,x);
		//Find kth element in xp
		double xk = gsl_vector_get(xp,k);
		//And replace with xk+delx:
		gsl_vector_set(xp,k,xk+delx);
		//Evaluate f(x+delx) there and store in fdx:
		f(xp,fdx);
		//Use fdx to calculate kth column of J:
		for(int i=0;i<dim;i++){
		//f(x+delx)
		double fdxi = gsl_vector_get(fdx,i);
		//f(x)
		double fxi  = gsl_vector_get(fx,i);
		//df/dx
		double fixk = (fdxi-fxi)/delx;
		//Store in J
		gsl_matrix_set(J,i,k,fixk);
		}
	}
	
	
	//Now we solve Jdelx = -f(x):
	
	
	//Since we have created the Jacobian matrix fdx is no longer used.
	//We now use it to store QR-decomposition in 
	//Copy content of J to Js, so that we dont destroy J
	gsl_matrix_memcpy(Js,J);
	gsl_vector_scale(fx,-1);//Correct the sign of f(x)
	//Decompose J<-QR and fdx<-tau
	gsl_linalg_QR_decomp(Js,fdx);
	//Solve system Js*delx = -fx and store xp <- delx
	gsl_linalg_QR_solve(Js,fdx,fx,xp);
	
	//Now we have delx stored in xp
	double lambda = 1;
	normfx = gsl_blas_dnrm2(fx);
	
	//Store vector x+delx in vector fx
	gsl_vector_memcpy(xpdx,x);
	gsl_blas_daxpy(lambda,xp,xpdx);
	//Store vector f(x+delx) in fdx:
	f(xpdx,fdx);
	double normfdelx = gsl_blas_dnrm2(fdx);	
	while(normfdelx>(1-lambda/2)*normfx && lambda>1./64){
	lambda /= 2;
	//Store vector x+delx in vector fx
	gsl_vector_memcpy(xpdx,x);
	gsl_blas_daxpy(lambda,xp,xpdx);
	//Store vector f(x+delx) in fdx:
	f(xpdx,fdx);
	normfdelx = gsl_blas_dnrm2(fdx);	
	}
	gsl_vector_memcpy(x,xpdx);

	//Recalculate f(x) and see if we have converged.
	f(x,fx);
	normfx = gsl_blas_dnrm2(fx);
	}

	gsl_matrix_free(J);
	gsl_vector_free(xp);
	gsl_vector_free(fx);
	gsl_vector_free(fdx);
	gsl_vector_free(xpdx);
	gsl_matrix_free(Js);

}
