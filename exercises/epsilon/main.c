#include<stdio.h>
#include<limits.h>
#include<float.h>
#include<math.h>

int equal(double a, double b, double tau, dobule epsilon){
	if (fabs(a-b) < tau)
		return 1;
	else if (fabs(a-b)/(fabs(a) + fabs(b)) < epsilon/2)
		return 1;
	else
		return 0;
}

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
	printf("\nNow we examine the machine epsilon:\n");

	double x=1; while(1+x!=1){x/=2;};x*=2;
	printf("machine eps (double and while)=%g\n",x);
	printf("machine eps =%g\n",DBL_EPSILON);

	float y=1; while(1+y!=1){y/=2;};y*=2;
	printf("machine eps (float and while)=%g\n",y);
	printf("machine eps =%g\n",FLT_EPSILON);
	
	long double z=1; while(1+z!=1){z/=2;};z*=2;
	printf("machine eps (long double and while)=%Lg\n",z);
	printf("machine eps =%Lg\n",LDBL_EPSILON);

	// Next part of the exercise:
	printf("\nNow we calculate some sums\n");
	int max = INT_MAX/2;
	float sum_up_float = 0;
	for(int i=1;i<=max;i++){
		sum_up_float+=1.0f/i;
	};
	printf("Limit of harmonic series (float) = %f\n",sum_up_float);

	float sum_down_float = 0;
	for(int i=0;i<max;i++){
		sum_down_float+=1.0f/(max-i);
	};
	printf("Limit of harmonic series (float) = %f\n",sum_down_float);
	// The last method is better, because float have few digits! Add like
	// numbers.
	
		
	double sum_up_double = 0;
	for(int i=1;i<=max;i++){
		sum_up_double+=1.0/i;
	};
	printf("Limit of harmonic series (double) = %g\n",sum_up_double);

	double sum_down_double = 0;
	for(int i=0;i<max;i++){
		sum_down_double+=1.0/(max-i);
	};
	printf("Limit of harmonic series (double) = %g\n",sum_down_double);
	// These two agree better due to higher precision.
	
	printf(\n"Now comes part 3. \n");

return 0;
}
