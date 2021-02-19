#include<stdio.h>
#include<math.h>
int main() {
	double x;
	int items;
	do{
		items=scanf("%lg",&x);
		if (items == EOF) break;
		printf("x=%g & sin(x) = %g & cos(x)=%g\n",x,sin(x),cos(x));
	}while(1);
return 0;
}
