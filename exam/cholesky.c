#include"cholesky_header.h"

//Printing functions have been written, mainly used for output
//(and while debugging)
void vector_print(gsl_vector * v){
	for(int i=0;i<v->size;i++){
		printf("%8g\n",gsl_vector_get(v,i));
	}
}

void matrix_print(gsl_matrix * M){
	int n=M->size1;
	int m=M->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		if(fabs(gsl_matrix_get(M,i,j))>1e-12){
		printf("%8g ",gsl_matrix_get(M,i,j));
		} else {
		printf("%8i ",0);
		}
		}
		printf("\n");
	}
	printf("\n");
}

void lineq_print(gsl_matrix * M,gsl_vector * v){
	int n=M->size1;
	int m=M->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
		printf("%8g ",gsl_matrix_get(M,i,j));
		}
		printf("x%i = %8g\n",i,gsl_vector_get(v,i));
	}
	printf("\n");
}

//This algorithm decomposes the real symmetric positive-definite matrix A into LL^T
void cholesky_decomp(gsl_matrix * A){
	//Determine dimension of problem
	int n=A->size1;
	int m=A->size2;
	//Assert that A is at least a square matrix? Just to avoid error...	
	assert(n==m);	
	//Sum over rows
	for(int i=0;i<n;i++){
		//Sum over entries in row left of diagonal	
		for(int j=0;j<i+1;j++){
			//Calculate the sum in the formula (41) in book
			double sum=0;
			for(int k=0;k<j;k++){
				double Lik = gsl_matrix_get(A,i,k);
				double Lij = gsl_matrix_get(A,j,k);
				sum += Lik*Lij;
			}
			//Special treatment for diagonal elements
			if(i==j){
				double Aii = gsl_matrix_get(A,i,i);
			       	gsl_matrix_set(A,j,j,sqrt(Aii-sum));
				for(int r=i+1;r<n;r++){
					gsl_matrix_set(A,i,r,0);
				}
			}
			//otherwise matrix elements are calculated like this:
			else{
				double Aij = gsl_matrix_get(A,i,j);
				double Ljj = gsl_matrix_get(A,j,j);
				double Lij = 1.0/Ljj*(Aij-sum);
				gsl_matrix_set(A,i,j,Lij);
			}
		}
	}
}

void rand_SPD(gsl_matrix * SPD){
	//This function creates a symmetric positive definite matrix 
	//by generating a diagonal matrix D with only positive entries
	//and a random matrix Q. Then the product Q'DQ is symmetric
	//positive definite.
	int n=SPD->size1;
	int m=SPD->size2;
	//Assert that matrix is square
	assert(n==m);
	
	//Allocate memory for D and Q
	gsl_matrix * D = gsl_matrix_alloc(n,m);
	gsl_matrix * Q = gsl_matrix_alloc(n,m);
	
	//Populate matrices randomly
	for(int i=0;i<n;i++){
		//Set eigenvalues to random positive number
		gsl_matrix_set(D,i,i,2*(double)rand()/RAND_MAX);
		for(int j=0;j<m;j++){
			//Populate random matrix Q
			gsl_matrix_set(Q,i,j,2*(double)rand()/RAND_MAX);
		}
	}
	//Calculate SPD = Q'DQ: First SPD=1*DQ
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,D,Q,0,SPD);
	//Then D = Q'SPD:
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Q,SPD,0,D);
	//Then copy D to SPD:
	gsl_matrix_memcpy(SPD,D);

	//And free temporary matrices:
	gsl_matrix_free(D);
	gsl_matrix_free(Q);
}

double cholesky_det(gsl_matrix * A){
	cholesky_decomp(A);
	//Calculate determinant of Cholesky-decomposition
	int n=A->size1;
	double determinant=1;
	for(int i=0;i<n;i++){
		//Just product of diagonal, since triangular matrix
		determinant*=gsl_matrix_get(A,i,i);
	}
	//Determinant of A is det(L)²:
	return determinant*determinant;
}

void cholesky_linsolve(gsl_matrix * A, gsl_vector * b, gsl_vector * x){
	int n=A->size1;
	//To solve a linear system Ax=b, use that A=LL':
	//To solve LL'x=b, note that L is triangular
	//Ly=b can then be solved by forward substitution
	//and L'x=y can then be solved by backwards substitution.
	
	//First, decompose A
	cholesky_decomp(A);
	//Then, solve Ly=b by forward substitution
	//(Note that in our notation it reads Ax=b)
	for(int i=0;i<n;i++){
		double sum = 0;
		for(int j=0;j<i;j++){
			sum+=gsl_vector_get(x,j)*gsl_matrix_get(A,i,j);
		}
		double xi = (gsl_vector_get(b,i)-sum)/gsl_matrix_get(A,i,i);
		gsl_vector_set(x,i,xi);
	}
	//Move y to vector named b
	gsl_vector_memcpy(b,x);
	//Solve L'x=y by backward substitution
	//(Note that it reads A'x=b in our notation)
	for(int i=n-1;i>=0;i--){
		double s = gsl_vector_get(b,i);
		for(int j=i+1;j<n;j++){
		//Note that indices on A are changed due to the transpose
			s-=gsl_matrix_get(A,j,i)*gsl_vector_get(x,j);
		}
		gsl_vector_set(x,i,s/gsl_matrix_get(A,i,i));
	}
	//Then x contains the solution to L'x=y
	//where y is solution to Ly=b, so that Ax = LL'x = b
}

void cholesky_inverse(gsl_matrix * A, gsl_matrix * B){
	//To find inverse, solve Ax_i = e_i, then {x1,x2,..} is inverse of A
	gsl_vector * e = gsl_vector_alloc(A->size1);
	gsl_vector * x = gsl_vector_alloc(A->size1);
	gsl_matrix * L = gsl_matrix_alloc(A->size1,A->size2);
	for(int i=0;i<A->size1;i++){
		gsl_vector_set(e,i,1);
		gsl_matrix_memcpy(L,A);
		cholesky_linsolve(L,e,x);
		gsl_matrix_set_col(B,i,x);
		gsl_vector_set_zero(e);
	}
	gsl_vector_free(e);
	gsl_vector_free(x);
	gsl_matrix_free(L);
}
