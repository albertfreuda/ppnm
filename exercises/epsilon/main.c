#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){
	int i=1; while(i+1>1) {i++;}
	printf("my max int (while) = %i\n",i-1);
	printf("compare with limits.h: %i\n",INT_MAX);
	i=1; do{
		i++;
		}while(i+1>1);
	printf("my max int (do while) = %i\n",i-1);
	for(int i=0;i+1>1;i++);
	printf("my max int (for) = %i\n",i-1);
	
	i=1;
	while(i-1<1) {i+=-1;}
	printf("my min int = %i\n",i+1);
	printf("compare with limits.h: %i\n",INT_MIN);
	i=1; do{
		i+=-1;
		}while(i-1<1);
	printf("my min int (do while) = %i\n",i+1);
	for(int i=0;i-1<1;i+=-1);
	printf("my min int (for) = %i\n",i+1);

	// Next part of the question
	printf("Now we examine the machine epsilon:\n");

	double x=1; while(1+x!=1){x/=2;};x*=2;
	printf("machine eps (double and while)=%g\n",x);
	printf("machine eps =%g\n",DBL_EPSILON);

	float y=1; while(1+y!=1){y/=2;};y*=2;
	printf("machine eps (float and while)=%f\n",y);
	printf("machine eps =%g\n",FLT_EPSILON);
	
	long double z=1; while(1+z!=1){z/=2;};z*=2;
	printf("machine eps (long double and while)=%Lg\n",z);
	printf("machine eps =%Lg\n",LDBL_EPSILON);
return 0;
}
