#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>


void print_matrix(gsl_matrix* A);

void vector_print(gsl_vector* x);

void GS_decomp(gsl_matrix* A,gsl_matrix* R){
	//Pull out the size of A:
	int n = A->size1, m = A->size2;
	//Make sure that A is long or quadratic:
	assert(n >= m);
	//Loop over columns of matrix A:
	for(int i=0;i<m;i++){
		//Calculate norm: Make room, add terms and squareroot it.
		double norm_ai=0;
		for(int j=0;j<n;j++){
	norm_ai += gsl_matrix_get(A,j,i)*gsl_matrix_get(A,j,i);
		}
		norm_ai = sqrt(norm_ai);
		//Assign to diagonals of R-matrix
		gsl_matrix_set(R,i,i,norm_ai);
		//Normalize a_i:
		for(int j=0;j<n;j++){
			double a_ji = gsl_matrix_get(A,j,i);
			if(norm_ai>1e-15){
			gsl_matrix_set(A,j,i,a_ji/norm_ai);}
			else {gsl_matrix_set(A,j,i,0);
			//	printf("Your matrix is singular!\n");
			}
		}
		//Orthogonalize all remaining vectors:
		for(int k=i+1;k<m;k++){
			double R_ik=0;
			//Calculate R's:
			for(int j=0;j<n;j++){
		R_ik += gsl_matrix_get(A,j,i)*gsl_matrix_get(A,j,k);
		gsl_matrix_set(R,i,k,R_ik);
			}
			//printf("R%i%i is %g\n",i+1,k+1,R_ik);
			//Calculate new orthogonal vectors:
			for(int j=0;j<n;j++){
		//printf("A%i%i = %g\n",j,k,gsl_matrix_get(A,j,k));
		//printf("A%i%i = %g\n",j,i,gsl_matrix_get(A,j,i));
		double a_jk = gsl_matrix_get(A,j,k)-gsl_matrix_get(A,j,i)*R_ik;
		gsl_matrix_set(A,j,k,a_jk);
		//printf("a%i%i = %g\n",j+1,k+1,a_jk);
			}
		}
	}
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	//First we change the system from QRx = b to Rx=Q^Tb
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
		//Then we perform backwards substitution:
	for(int i=x->size-1;i>=0;i--){
		double s = gsl_vector_get(x,i);
		for(int k=i+1;k<x->size;k++){
			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
		}
		if(fabs(s/gsl_matrix_get(R,i,i))>1e-12){
			gsl_vector_set(x,i,s/gsl_matrix_get(R,i,i));
		} else {
			gsl_vector_set(x,i,0);	
		}
		
	}
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	//To find inverse, we solve Axi=ei for all i. Then {x1,x2,x3,...}
	//is the inverse of A.
	// Make room for unit vector:
	gsl_vector* e = gsl_vector_alloc(R->size2);
	gsl_vector* x = gsl_vector_alloc(R->size2);
	for(int i=0;i<R->size2;i++){
		//Create the proper unit vector:
		gsl_vector_set(e,i,1);
		//Solve the system and save to x:
		GS_solve(Q,R,e,x);
		//Save x to column of matrix B:
		for(int k=0;k<R->size2;k++){
			gsl_matrix_set(B,k,i,gsl_vector_get(x,k));
		}
		gsl_vector_set(e,i,0);		
	}	
}






