#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>

void vector_print(const char* s,gsl_vector * v);

void reflection(gsl_vector* hi, gsl_vector* ce){
	gsl_vector_scale(hi,-1);
	gsl_blas_daxpy(2,ce,hi);
}

void expansion(gsl_vector * hi, gsl_vector* ce){
	gsl_vector_scale(hi,-2);
	gsl_blas_daxpy(2,ce,hi);
}

void contraction(gsl_vector * hi,gsl_vector * ce){
	gsl_vector_scale(hi,0.5);
	gsl_blas_daxpy(0.5,ce,hi);
}
	
void reduction(gsl_matrix * simplex,
	       gsl_vector * lo, 
	       gsl_vector * x,
	       int kmin,
	       int n){
	
	for(int i=0;i<n+1;i++){
		if(i!=kmin){
		gsl_matrix_get_col(x,simplex,i);
		gsl_vector_add(x,lo);
		gsl_vector_scale(x,0.5);	
		gsl_matrix_set_col(simplex,i,x);
		}
	}
}

double size_of_simplex(gsl_matrix * S){
	int n=S->size1;

	double s = 0,dist;
for(int k=0;k<n;k++){
	for(int i=k+1;i<n+1;i++){
		dist = 0;
		for(int l=0;l<n;l++){
		double pk = gsl_matrix_get(S,l,k);
		double pi = gsl_matrix_get(S,l,i);
		dist+=(pk-pi)*(pk-pi);
		}
		if(dist>s*s) s=sqrt(dist);
	}
}
return s;
}

void downhill_simplex(double f(gsl_vector * x),
		gsl_matrix * simplex,
		double tolerance){
	//Name dimension of problem n:
	int n = simplex->size1;
	
	//Find highest,lowest and centroid of simplex
	gsl_vector* x  = gsl_vector_alloc(n);//Make room for some point
	gsl_vector* hi = gsl_vector_alloc(n);//Make room for highest
	gsl_vector* lo = gsl_vector_alloc(n);//Make room for lowest
	gsl_vector* ct = gsl_vector_alloc(n);//Make room for centroid
	gsl_vector* reflected  = gsl_vector_alloc(n);//Make room for centroid
	gsl_vector* expanded   = gsl_vector_alloc(n);//Make room for centroid
	gsl_vector* contracted = gsl_vector_alloc(n);//Make room for centroid
	gsl_vector* reduced    = gsl_vector_alloc(n);//Make room for centroid
	
	//Determine f-value of first vertex:
	gsl_matrix_get_col(hi,simplex,0);//Copy 0th column into hi
	       	gsl_vector_memcpy(lo,hi);
	double fmax = f(hi),newval,fmin = f(lo);	
	
	int kmax = 0,kmin=0;
//Run through all vectors in simplex to find hi and lo:
for(int k=1;k<n+1;k++){
	gsl_matrix_get_col(x,simplex,k);//Copy kth column into x
	newval = f(x);//Calculate f of next vertex
       	if(newval>fmax){//If larger, this is max
	       	gsl_vector_memcpy(hi,x);
		fmax = newval;kmax=k;
	}
	if(newval<fmin){
		gsl_vector_memcpy(lo,x);
		fmin = newval;kmin = k;
	}
}
//Find centroid: First add all vertices except max
for(int k=0;k<n+1;k++){
	if(k!=kmax){
	gsl_matrix_get_col(x,simplex,k);//Copy kth column into x
	gsl_vector_add(ct,x);
	}
}
//Then scale to find centroid:
gsl_vector_scale(ct,1./n);

double size = size_of_simplex(simplex);

//vector_print("HI:",hi);
//vector_print("LO:",lo);
//vector_print("CT:",ct);

//Now start moving:
while(size>tolerance){
//Perform reflection and save reflected point in reflected	
gsl_vector_memcpy(reflected,hi);
reflection(reflected,ct);
if(f(reflected)<f(lo)){
	//Try expansion:
	gsl_vector_memcpy(expanded,hi);
	expansion(expanded,ct);
	if(f(expanded)<f(reflected)){
		//If expansion was succesful:
		//Replace highest vertex in simplex:
		gsl_matrix_set_col(simplex,kmax,expanded);
		//printf("expanded\n");
	//Otherwise, switch hi for reflected:
	} else gsl_matrix_set_col(simplex,kmax,reflected);
		//printf("reflected (smaller than lo)\n");
	
}
else{
	if(f(reflected)<f(hi)){
		gsl_matrix_set_col(simplex,kmax,reflected);
		//printf("reflected\n");
	}
	else{
		//Try contraction
		gsl_vector_memcpy(contracted,hi);
		contraction(contracted,ct);
		if(f(contracted)<f(hi)){
			gsl_matrix_set_col(simplex,kmax,contracted);
			//printf("contracted\n");
		} 
		else {reduction(simplex,lo,x,kmin,n);
		       	//printf("reduced\n");
		}
	}
}
	
//Determine f-value of first vertex:
gsl_matrix_get_col(hi,simplex,0);//Copy 0th column into hi
	gsl_vector_memcpy(lo,hi);
fmax = f(hi);
fmin = f(lo);	

kmax = 0,kmin=0;
//Run through all vectors in simplex to find hi and lo:
for(int k=1;k<n+1;k++){
	gsl_matrix_get_col(x,simplex,k);//Copy kth column into x
	newval = f(x);//Calculate f of next vertex
       	if(newval>fmax){//If larger, this is max
	       	gsl_vector_memcpy(hi,x);
		fmax = newval;kmax=k;
	}
	if(newval<fmin){
		gsl_vector_memcpy(lo,x);
		fmin = newval;kmin = k;
	}
}
//Find centroid: First add all vertices except max
gsl_vector_set_zero(ct);
for(int k=0;k<n+1;k++){
	if(k!=kmax){
	gsl_matrix_get_col(x,simplex,k);//Copy kth column into x
	gsl_vector_add(ct,x);
	}
}
//Then scale to find centroid:
gsl_vector_scale(ct,1./n);

size = size_of_simplex(simplex);
//printf("Size = %g\n",size);
//vector_print("HI:",hi);
//vector_print("LO:",lo);
//vector_print("CT:",ct);
}

	gsl_vector_free( x);
	gsl_vector_free(hi);
	gsl_vector_free(lo);
	gsl_vector_free(ct);
	gsl_vector_free( reflected);
	gsl_vector_free(  expanded);
	gsl_vector_free(contracted);
	gsl_vector_free(   reduced);
}

