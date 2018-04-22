#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Runge_Kutta.h"
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
	// part A and B

	int n=2, k, max=1e6; double b=7, h=b/20, acc=1e-26, eps=1e-2;
	void system1(int n,double x,gsl_vector* y,gsl_vector* dydx){
		double y0=gsl_vector_get(y,0);
		double y1=gsl_vector_get(y,1);
		double f0=y1*x;
		double f1=-y0;
		gsl_vector_set(dydx,0,f0);
		gsl_vector_set(dydx,1,f1);		
	}
	void orbital_equation(int n, double t,gsl_matrix* y, gsl_matrix* dydt) {
	
	(void)(t); 
	double eps = 0.01;
	double y0=gsl_vector_get(y,0);
	double y1=gsl_vector_get(y,1);
	gsl_vector_set(dydt,0,y1);
	gsl_vector_set(dydt,1,1-y0+eps*y0*y0);
	}


	gsl_vector* xlist=gsl_vector_alloc(max+1);
	gsl_matrix* ylist=gsl_matrix_alloc(max+1,n);
	gsl_vector* y=gsl_vector_alloc(n);

	printf("A. rk21 midpoint euler stepper \n\n");
	// initialize
	gsl_vector_set(xlist,0,0.0);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.0);
	
	printf("Testing on the system \n");
	printf("dydx_0= y1*x\n");
	printf("dydx_1= -y0\n");
	printf("eps=%g, acc=%g, max=%i, b=%g, hstart=%g\n",eps,acc,max,b,h);
	printf("Startpoint set to:\n");
	printf("x=%g",gsl_vector_get(xlist,0));
	printf("y set to:\n");
	gsl_matrix_get_row(y,ylist,0);
	printv(y);
	k=ode_driver(system1,n,xlist,ylist,b,h,acc,eps,max);
	printf("And the solution at y(%g):\n",b);
	gsl_matrix_get_row(y,ylist,k);
	printv(y);
	printf("Number of steps taken: %i\n\n",k);


	h=1e-3; acc=1e-6; eps=1e-6; b=M_PI*20; max=1e6;
	gsl_vector_set_zero(xlist);
	gsl_matrix_set_zero(ylist);
	gsl_vector_set(xlist,0,0.0);
	gsl_matrix_set(ylist,0,0,1.0);
	gsl_matrix_set(ylist,0,1,0.5);
	printf("Testing on the system from the orbital equation exercise in PP\n");
	printf("with eps=0.01 (relativistic correction)\n");
	printf("dydx_0= y1\n");
	printf("dydx_1= 1-y0+eps*y0*y0\n");
	printf("eps=%g, acc=%g, max=%i, b=%g, hstart=%g\n",eps,acc,max,b,h);
	printf("Startpoint set to:\n");
	printf("x=%g",gsl_vector_get(xlist,0));
	printf("y set to:\n");
	gsl_matrix_get_row(y,ylist,0);
	printv(y);
	k=ode_driver(orbital_equation,n,xlist,ylist,b,h,acc,eps,max);
	printf("And the solution at y(%g):\n",b);
	gsl_matrix_get_row(y,ylist,k);
	printv(y);
	printf("Number of steps taken: %i\n",k);
	printf("The value I got for y0 in the original exercise was 0.700198\n");
	printf("So this is very close. That was of course using gsl library rk8pd\n");
	printf("algorithm which probably is more acurate than the simple implementation here\n");











	gsl_vector_free(xlist);
	gsl_matrix_free(ylist);
	gsl_vector_free(y);

	return EXIT_SUCCESS;
}