#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Runge_Kutta.h"
#include "Integrator.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.6f" //format of print "7 width, 3 digits after comma" 


// inspired by Dmitri Fedorovs print implementation
void printm(gsl_matrix* A) {
	// iterating over all rows and columns 
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++) {
			printf(FMT,gsl_matrix_get(A,i,j));
		}
		printf("\n");}
}

void printv(gsl_vector *A){
	for(int i=0;i<A->size;i++){
		printf(FMT,gsl_vector_get(A,i));
		printf("\n");
	}
}

int main(void) {
	// part A 

	double a=0,b=1,acc=1e-3,eps=1e-6; int calls=0;

	double f1(double x) {calls++; return sqrt(x);}
	double f2(double x) {calls++; return 1.0/sqrt(x);}
	double f3(double x) {calls++; return log(x)/sqrt(x);}
	double f4(double x) {calls++; return 4*sqrt(1-(1-x)*(1-x));}
	

	double result1 = adapt(f1,a,b,acc,eps); 



	printf("A. Recursive adaptive integration\n\n");
	
	printf("Integrating sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=2/3, Calculated=%g\n",result1);
	printf("Number of function calls: %i\n\n",calls); calls=0;
	double result2 = adapt(f2,a,b,acc,eps); 
	printf("Integrating 1/sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=2, Calculated=%g\n",result2);
	printf("Number of function calls: %i\n\n",calls);calls=0;
	double result3 = adapt(f3,a,b,acc,eps); 
	printf("Integrating ln(x)/sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=-4, Calculated=%g\n",result3);
	printf("Number of function calls: %i\n\n",calls);calls=0;
	eps=1e-18; acc=1e-18; 
	double result4 = adapt(f4,a,b,acc,eps);
	printf("Integrating 4*sqrt(1-(1-x)^2) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=pi, Calculated=%7.20f\n",result4);
	printf("Number of function calls: %i\n\n",calls);



	return EXIT_SUCCESS;
}