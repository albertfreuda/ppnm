#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
	double x = tgamma(5);
	printf("Gamma function of 5 is: %g\n",x);
	
	double y = j1(0.5);
	printf("Bessel function of 0.5 is: %g\n",y);

	complex double z = csqrt(-2);
	printf("Sqrt(-2) is: %g + %g i\n",creal(z),cimag(z));

	complex double a = cexp(I*M_PI);
	printf("e^(i*pi) is: %g + %g i\n",creal(a),cimag(a));

	complex double b = cexp(I);
	printf("e^(i) is: %g + %g i\n",creal(b),cimag(b));

	complex double c = cpow(I,M_E);
	printf("i^e is: %g + %g i\n",creal(c),cimag(c));
	
	complex double d = cpow(I,I);
	printf("i^i is: %g + %g i\n",creal(d),cimag(d));

	double e_double = 1.0/9;
	float e_float = 1.f/9;
	long double e_long_double = 1.L/9;


	printf("float is: %.25g\n",e_float);
	printf("double is: %.25lg\n",e_double);
	printf("long_double is: %.25Lg\n",e_long_double);
	return 0;
}
