#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void print_matrix(gsl_matrix* A);
void print_vector(gsl_vector* x);

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta), s = sin(theta);
	for(int i=0;i<A->size1;i++){
		double a_ip = c*gsl_matrix_get(A,i,p) - s*gsl_matrix_get(A,i,q);
		double a_iq = s*gsl_matrix_get(A,i,p) + c*gsl_matrix_get(A,i,q);
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
	gsl_matrix_set_identity(V);
	int changed;
	do{
		changed = 0;
		for(int p=0;p<A->size1;p++){
	for(int q=p+1;q<A->size1;q++){
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
		//Check if the new element is different from the old element:
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

void optim_jacobi_diag(gsl_matrix* A,gsl_vector* e, gsl_matrix* V){
	//Two (quadratic) matrices, decomposes A into VDV'
	//A is overwritten; A<-D so the matrix A is returned as diagonal
	//V must be the identity, so we make it:
	gsl_matrix_set_identity(V);
	int changed,n=A->size1;
	for(int i=0;i<n;i++) gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
	do{
		changed = 0;
		for(int p=0;p<n;p++){
	for(int q=p+1;q<n;q++){
		//Extract relevant elements in A
		double A_pq = gsl_matrix_get(A,p,q);
		double A_pp = gsl_vector_get(e,p);
		double A_qq = gsl_vector_get(e,q);
		//Calculate Jacobi rotation angle:
		double theta = 0.5*atan2(2*A_pq,A_qq-A_pp);
		double c = cos(theta),s = sin(theta);//To shorten notation
		//Calculate new diagonal elements. 
		double new_App = c*c*A_pp - 2*s*c*A_pq+s*s*A_qq;
		double new_Aqq = s*s*A_pp + 2*s*c*A_pq+c*c*A_qq;
		//Check if the new element is different from the old element:
		if(new_App!=A_pp || new_Aqq!=A_qq){
			//If the diagonals actually changed, we haven't converged
			//Perform rotation again:
	changed=1; //Update program status
	gsl_vector_set(e,p,new_App);//Save updated eigenvalues
	gsl_vector_set(e,q,new_Aqq);
	gsl_matrix_set(A,p,q,0.0);//Set off-diagonal to zero (by construction)
	//A <- J'AJ
	for(int i=0;i<p;i++){//Update down to p-diag
		//We don't need to touch diagonal, so i<p
		double aip = gsl_matrix_get(A,i,p);
		double aiq = gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,c*aip-s*aiq);
		gsl_matrix_set(A,i,q,s*aip+c*aiq);}
	for(int i=p+1;i<q;i++){//Update inside Jacobi square
		// We don't need to change Apq, so i=p+1
		double api = gsl_matrix_get(A,p,i);//Api = Aip
		double aiq = gsl_matrix_get(A,i,q);//Above diagonal, so it's fine
		gsl_matrix_set(A,p,i,c*api-s*aiq);
		gsl_matrix_set(A,i,q,c*aiq+s*api);}
	for(int i=0;i<p;i++){//Update right and out
		double api = gsl_matrix_get(A,p,i);	
		double aqi = gsl_matrix_get(A,q,i);
		gsl_matrix_set(A,p,i,c*api-s*aqi);
		gsl_matrix_set(A,q,i,c*aqi+s*api);}	
	for(int i=0;i<n;i++){//Update V-matrix
		double vip = gsl_matrix_get(V,i,p);
		double viq = gsl_matrix_get(V,i,q);
		gsl_matrix_set(V,i,p,c*vip-s*viq);
		gsl_matrix_set(V,i,q,c*viq+s*vip);}	
		}
	}
		}
		}while(changed!=0);
}
