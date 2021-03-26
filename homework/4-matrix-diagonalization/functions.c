#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta), s = sin(theta);
	for(int i=0;i<A->size1;i++){
		double a_ip = c*gsl_matrix_get(A,i,p) - s*gsl_matrix_get(A,i,q);
		double a_iq =-s*gsl_matrix_get(A,i,p) + c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,a_ip);
		gsl_matrix_set(A,i,q,a_iq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta), s = sin(theta);
	for(int i=0;i<A->size1;i++){
		double a_pi = c*gsl_matrix_get(A,p,i) + s*gsl_matrix_get(A,q,i);
		double a_qi =-s*gsl_matrix_get(A,p,i) + c*gsl_matrix_get(A,q,i);
		gsl_matrix_set(A,p,i,a_pi);
		gsl_matrix_set(A,q,i,a_qi);
	}
}

void jacobi_diag(gsl_matrix* A,gsl_matrix* V){
	//Two (quadratic) matrices, decomposes A into VDV'
	//A is overwritten; A<-D so the matrix A is returned as diagonal
	//V must be the identity, so we make it:
	for(int i=0;i<V->size1;i++){
		gsl_matrix_set(V,i,i,1);
	}
	int changed;
	do{
		changed = 0;
		for(int p=0;p<A->size1-1;p++){
	for(int q=p+1;q<A->size2;q++){
		//Extract relevant elements in A
		double A_pq = gsl_matrix_get(A,p,q);
		double A_pp = gsl_matrix_get(A,p,p);
		double A_qq = gsl_matrix_get(A,q,q);
		//Calculate Jacobi rotation angle:
		double theta = 0.5*atan2(2*A_pq,A_qq-A_pp);
		double c = cos(theta),s = sin(theta);//To shorten notation
		//Calculate new diagonal elements. 
		double new_App = c*c*A_pp - 2*s*c*A_pq+s*s*A_qq;
		double new_Aqq = s*s*A_pp + 2*s*c*A_pq+c*c*A_qq;
		printf("%i%i and theta = %g\n",p,q,theta);
		//Check if the new element is different from the old element:
		//printf("pq = %i%i: Theta = %g, and we have Apq = %g\n",p,q,theta,A_pq);
		if(new_App!=A_pp || new_Aqq!=A_qq){
			//If the diagonals actually changed, we haven't converged
			//Perform rotation again:
			changed=1; //Update program status
			timesJ(A,p,q, theta);//Update A<-J'AJ
			Jtimes(A,p,q,-theta);//Update A (note J'(t) = J(-t))
			timesJ(V,p,q, theta);//Update V<-VJ
		}
	}
		}
		}while(changed!=0);
}
