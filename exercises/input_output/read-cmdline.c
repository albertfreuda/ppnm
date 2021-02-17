#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// The purpose of this program is to read a set of numbers and print them,
// along with their sine and cosines to the stdout.

int main(int argc, char** argv){
	printf("Thanks for running my %s program!\n",argv[0]);
	for(int i=1;i<argc;i++){
		double x = atof(argv[i]);
		printf("x = %g, sin(x) = %g, cos(x) = %g\n",x,sin(x),cos(x));
	}
	return 0;
}
