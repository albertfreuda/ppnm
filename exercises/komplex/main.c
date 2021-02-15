#include"komplex.h"
#include<stdio.h>
#define TINY 1e-6

int main(){
	komplex a = {1,2}, b = {3,4};
	
	printf("Here we test the functions in komplex.\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should equal ",R);
	komplex_print("when in fact they equal",r);

	return 0;
}

