----------------------------------INTRO--------------------------------------
Hello, and welcome to my exam project. My student number ends with 95, so

95 mod 22 = 7

meaning that I have implemented the Cholesky-Banachiewicz algorithm, 
for performing a Cholesky decomposition of a symmetric positive definite 
real matrix A into LL', where L is a lower triangular matrix. The algorithm 
does so 'in place' to save memory. The signature is:

void cholesky_decomp(gsl_matrix * A)

The matrix objects from the GNU Scientific Library are used as containers. 
The original matrix is destroyed, and replaced with the lower triangular 
matrix L, such that LL' = A.

----------------------------------TESTING-------------------------------------
To test the algorithm, a random symmetric positive definite matrix must be 
produced. This is done by creating a diagonal matrix D with positive diagonals 
and a random matrix Q. Then:

A = QDQ'

is a symmetric positive definite matrix. This function has signature

void rand_SPD(gsl_matrix * A)

and creates the symmetric positive definite matrix in A.

--------------------------------DETERMINANT-----------------------------------
A determinant function has also been implemented. 

double cholesky_det(gsl_matrix * A)

This function takes a real symmetric positive definite matrix, performs the
cholesky-decomposition via cholesky_decomp, and calculates the determinant as

det(A) = det(LL') = det(L)det(L')=det(L)²

Since L is triangular, the determinant is simply the product of diagonal entries:

det(L) = L(1,1)*L(2,2)*...L(n,n)

Since the cholesky_det function performs the decomposition, it destroys the
matrix when running. 

-----------------------------LINEAR-EQUATION-SOLVER--------------------------
A linear equation solver has also been implemented. It solves a linear system 

Ax = b

when A is real, symmetric and positive definite. It does so by first performing
the Cholesky decomposition using cholesky_decomp (the matrix is destroyed):

A = LL'

and then solving LL'x=b in two steps: The system Ly=b can be solved by forward
substitution, because L is lower triangular. Then, the system L'x=y can be solved
by backward substition (because L' is upper triangular). In total

b = Ly = L(L'x) = LL'x = Ax

---------------------------------MATRIX-INVERSE--------------------------------
To calculate the matrix inverse, I solve n systems of linear equations:

Ax_i = e_i     (1)

Then, from the definition of matrix multiplication, it follows that the matrix

B = {x1,x2,x3,...}

must obey

AB = I

I solve the linear system in (1) by using the cholesky_linsolve function. 
