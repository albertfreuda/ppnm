#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void print_vector(gsl_vector* v){
	int n = v->size;
	for(int i=0;i<n;i++){
		printf("%g\n",gsl_vector_get(v,i));
	}
}

void print_matrix(gsl_matrix* A){
	int n = A->size1;
	int m = A->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		printf("%g ",gsl_matrix_get(A,i,j));
		}
		printf("\n");	
	}
}

//We create the stepper-function:
void rkstep12(
	void f(double t,gsl_vector* y,gsl_vector* dydt),//dydt=f(t,y)
	double t,                                       //current value	t
	gsl_vector* yt,                                 //current value y(t)
	double h,                                       //step size h
	gsl_vector* yh,					//out: y(t+h)
	gsl_vector* err){				//out: error estimate dely

	int n=yt->size;
	gsl_vector* k0 = gsl_vector_alloc(yt->size);
	gsl_vector* k  = gsl_vector_alloc(yt->size);
	gsl_vector* ymid = gsl_vector_alloc(yt->size);

	//Take step and evaluate. Store in k0
	f(t,yt,k0);
	
	//Use k0 to take half-step:
	for(int i=0;i<n;i++){
	gsl_vector_set(ymid,i,gsl_vector_get(yt,i)+0.5*h*gsl_vector_get(k0,i));
	}
	
	f(t+0.5*h,ymid,k);//Use estimate to evaluate at midpoint. Save to k
	for(int i=0;i<n;i++) gsl_vector_set(yh,i,gsl_vector_get(yt,i)+h*gsl_vector_get(k,i));
	for(int i=0;i<n;i++) gsl_vector_set(err,i,(gsl_vector_get(k0,i)-gsl_vector_get(k,i))/2);
	gsl_vector_free(k0);
	gsl_vector_free(k);
	gsl_vector_free(ymid);
}

void driver(void f(double t, gsl_vector* y,gsl_vector* dydt),
		double a,//t0
		gsl_vector* ya, //y(a)
		double b,//end point of integration
		gsl_vector* yb, //y(b) to be calc
		gsl_vector* err, //err memory location
		double h,//initial step size
		double acc,//abs acc. goal
		double eps){//relative acc. goal
		
	//Calculate y(b)
	double t = a;//Starting point
	printf("%g ",t); //Print first t and initial condition y(a):
	for(int i=0;i<ya->size;i++){
		printf("%g ",gsl_vector_get(ya,i));
	}	
	printf("\n");
	while(t<b){//Check that we are still in the interval
		if(t+h>b) h = b-t;//Don't step too far if we reach end point
		//Take step:
		rkstep12(f,t,ya,h,yb,err);
		//New y(t+h) is stored in yb
		//Check if tolerance is satisfied:
		double norm_yb = gsl_blas_dnrm2(yb);
		double norm_err = gsl_blas_dnrm2(err);
		double tol = (eps*norm_yb+acc)*sqrt(h/(b-a));
		//Check if norm_err<tol
		if(norm_err<tol){ //If true accept step
			t = t+h; //Progress time-step
			gsl_vector_memcpy(ya,yb);//Save new step to current step
			//Print out new point:
			printf("%g ",t);
			for(int i=0;i<ya->size;i++){
				printf("%g ",gsl_vector_get(ya,i));
			}	
			printf("\n");
		}
		//Calculate new step-size h
		h *= 0.95*pow(tol/norm_err,0.25);
	}
}
