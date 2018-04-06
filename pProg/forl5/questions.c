#include <stdio.h>
#include <stdlib.h>
int main() {
	struct vector {double x,y,z;};
	struct vector v = {1,2,3};
	struct vector x = {1,2,3};
	struct vector z = {1,2,3};
	struct vector u = {.x=1,.y=2,.z=3};
	struct vector w = {.z=3,.y=2,.x=1};
	
	typedef struct vector vector;
	vector o;

  	double* a=foo(); a[2]=1;
  	double* b=foo(); b[2]=2;
  	printf("a[2] = %g\n",a[2]);
	return EXIT_SUCCESS;
}

double* foo(){ double a[]={0,0,0}; return &a; }