Hello, and welcome to my exam project. I have implemented the
Cholesky-Banachiewicz algorithm, for performing a Cholesky decomposition 
of a symmetric positive definite real matrix A into LL', where L is a 
lower triangular matrix. The algorithm does so 'in place' to save memory. 
The signature is:

void cholesky_decomp(gsl_matrix * A)

The matrix objects from the GNU Scientific Library are used as containers. 
The original matrix is destroyed.

To test the algorithm, a random symmetric positive definite matrix must be 
produced. This is done by creating a diagonal matrix D with positive diagonals 
and a random matrix Q. Then:

A = QDQ'

is a symmetric positive definite matrix. This function has signature

void rand_SPD(gsl_matrix * A)

A determinant function has also been implemented. 

double cholesky_det(gsl_matrix * A)

This function takes a real symmetric positive definite matrix, performs the
cholesky-decomposition vi cholesky_decomp, and calculates the determinant as

det(A) = det(LL') = det(L)det(L')=det(L)²

Since L is triangular, the determinant is simply the product of diagonal entries.

det(L) = L(1,1)*L(2,2)*...L(n,n) 
