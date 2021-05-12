#include"minimization.h"
#define delx  sqrt(DBL_EPSILON)

void vector_print(const char* s,gsl_vector * v){
	printf("%s ( ",s);
	int n=v->size;
	for(int i=0;i<n;i++){
		printf("%.3g ",gsl_vector_get(v,i));
	}
	printf(")\n");
}
void matrix_print(const char* s,gsl_matrix * M){
	printf("%s\n",s);
	int n=M->size1,m=M->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		printf("%.3g ",gsl_matrix_get(M,i,j));
		}
		printf("\n");
	}
}

//Gradient function seems to work perfectly fine. Tested
void ngradient(double F(gsl_vector * x),
		gsl_vector * x,
		gsl_vector * grad)
		{
	int dim = x->size;
	for(int i = 0;i<dim;i++){
		//Store function value:
		double f_of_x = F(x);
		//Store i'th x-value (we'll change it in a bit)
		double dx,xi = gsl_vector_get(x,i);
		//Since F might not be normalized
		//we correct the step-size accordingly:
		if(fabs(xi)<delx) dx = delx;
		else dx = fabs(xi)*delx;	
		//Move delta x in ith direction:
		gsl_vector_set(x,i,xi+dx);
		//Evaluate at new x:
		double f_of_xpdx = F(x);
		//Calculate gradient:
		double dFdxi = (f_of_xpdx - f_of_x) / dx;
		//Store result:
		gsl_vector_set(grad,i,dFdxi);
		//Go back to be ready for next dimension:
		gsl_vector_set(x,i,xi);
	}	
}

int qnewton(double F(gsl_vector * x), //function to minimize
		gsl_vector * x,        //initial guess
		double tol)            //tolerance
		{
int dim = x->size,nsteps=0;
gsl_matrix * B = gsl_matrix_alloc(dim,dim);
gsl_vector * deltax   = gsl_vector_alloc(dim);
gsl_vector * gradient = gsl_vector_alloc(dim);
gsl_vector * z        = gsl_vector_alloc(dim);
gsl_vector * u        = gsl_vector_alloc(dim);

//Set the inverse Hessian to identity as a start guess:
gsl_matrix_set_identity(B);

while(nsteps<1000){
	nsteps++;
	//First we calculate delta x:
ngradient(F,x,gradient);//Calculate gradient at x
gsl_blas_dgemv(CblasNoTrans,-1,B,gradient,0,deltax);//dx = -B*gradF(x)

//------------------------------------------------------------
//matrix_print("The inverse of the hessian is:",B);
//vector_print("Delta x is calculated to be:",deltax);

	//Now we perform linesearch
//Set lambda=1 to begin:
double lambda = 1,lambda_min = DBL_EPSILON;
gsl_vector_memcpy(z,x);//Store copy of point
gsl_vector_add(z,deltax);//Add step to copy
double fz = F(z),fx = F(x);//calculate function at old and new point
double sTg;
gsl_blas_ddot(deltax,gradient,&sTg);//Caculate important sTg

//-----------------------------------------------------------
//printf("Before linesearch fx=%g and fz=%g\n",fx,fz);

while(fx+0.01*sTg<=fz){
	//If we reach lambda_min, reset hessian to identity
	if(lambda<lambda_min) {
	//	printf("Lambda=%g, and hessian is reset.\n",lambda);
		gsl_matrix_set_identity(B); break;
		}
	lambda *= 0.5;//Smaller step
	gsl_vector_scale(deltax,0.5);//Smaller step
	gsl_vector_memcpy(z,x);//Step back
	gsl_vector_add(z,deltax);//Add smaller step to copy
	fz = F(z);//calculate function at new point
	gsl_blas_ddot(deltax,gradient,&sTg);//Caculate importan sTg
}
//----------------------------------------------------------
//vector_print("After linesearch our step is:",deltax);
//printf("and our fz=%g\n",fz);

	//Update our point:
gsl_vector_memcpy(x,z);
//vector_print("Best estimate of minimum is:",x);
	
	//Calculate gradient at new x<-x+dx
ngradient(F,x,z);//store new gradient in z
	
	//See if accuracy is met:
if(gsl_blas_dnrm2(z)<tol){
//	printf("Hooray, we have met the accuracy goal at fx = %g !\n",fz);
	break;
}

	//Now we calculate update to inverse Hessian:
gsl_vector_sub(z,gradient);//z<-grad(x+dx)-grad(x) (z=y in note)
gsl_vector_memcpy(u,deltax);//copy deltax to u
gsl_blas_dgemv(CblasNoTrans,-1,B,z,1,u);//u <- -B*y+deltax
double sTy,uTy;
gsl_blas_ddot(deltax,z,&sTy);//Calculate scaling of u
//If scaling is not outrageous, update inverse hessian. Else, leave it
if(fabs(sTy)>1e-6){
	gsl_blas_ddot(u,z,&uTy);//Calculate scaling of u
	double gamma=uTy/2/sTy;
	gsl_blas_daxpy(-gamma,deltax,u);
	gsl_blas_dger(1./sTy,u,deltax,B);//B <- B + u*sT/(sT*y)
	gsl_blas_dger(1./sTy,deltax,u,B);//B <- B + u*sT/(sT*y)
}
//---------------------------------------------------------
//else printf("Ouch, the update is quite big - so we don't apply it\n");
}


gsl_matrix_free(B);
gsl_vector_free(deltax);
gsl_vector_free(gradient);
gsl_vector_free(z);
gsl_vector_free(u);
return nsteps;
}
