#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main (int argc, char *argv[]){
	char* input_name = argv[1];
	char* output_name = argv[2];
	
	FILE* my_input_stream = fopen(input_name,"r");
	FILE* my_output_stream = fopen(output_name,"w");
	double x;
	int items;
	do{
		items=fscanf(my_input_stream,"%lg",&x);
		fprintf(my_output_stream,"x=%g, sin(x)=%g, cos(x)=%g\n",x,sin(x),cos(x));
	}while(items!=EOF);
	
	fclose(my_input_stream);
	fclose(my_output_stream);
	return 0;
}

