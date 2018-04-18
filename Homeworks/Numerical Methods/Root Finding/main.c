#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Newton.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.8f" //format of print "7 width, 3 digits after comma" 


// inspired by Dmitri Fedorovs print implementation
void printm(gsl_matrix *A){
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

int main() {
	// part A
	int n=2;

	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);

	void system_f(gsl_vector* x, gsl_vector* fx){
	double A = 10000.0;
	double x0 = gsl_vector_get(x,1);
	double y =	gsl_vector_get(x,0.1);
	double fx0 = A*x0*y -1; 
	double fx1 = exp(-x0) + exp(-y) -(1.0 + 1.0/A);
	gsl_vector_set(fx,0,fx0);
	gsl_vector_set(fx,1,fx1);
	}

	void Rosenbrock_grad_f(gsl_vector* x, gsl_vector* fx){
	double x0 = gsl_vector_get(x,0);
	double y =	gsl_vector_get(x,1);
	double fx0 = 2*x0-2+400*(x0*x0*x0-y*x0); 
	double fx1 = 200*(y-x0*x0);
	gsl_vector_set(fx,0,fx0);
	gsl_vector_set(fx,1,fx1);
	}

	void Himmelblau_grad_f(gsl_vector* x, gsl_vector* fx){
	double x0 = gsl_vector_get(x,0);
	double y =	gsl_vector_get(x,1);
	double fx0 = 4*x0*(x0*x0+y-11)+2*(x0+y*y-7); 
	double fx1 = 2*(x0*x0+y-11)+4*y*(x0+y*y-7);
	gsl_vector_set(fx,0,fx0);
	gsl_vector_set(fx,1,fx1);
	}
	
	
	// system
	double dx=1e-9, epsilon=1e-3;
	printf("Testing if the implementation works on the system of equations provided.\n");
	
	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton_num(system_f,x,dx,epsilon);	
	
	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is (should be close to zero-vector)\n");
	system_f(x,fx);
	printv(fx);

	// Rosenbrock
	dx=1e-8, epsilon=1e-4;
	printf("Testing if the implementation works on Rosenbrock.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton_num(Rosenbrock_grad_f,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Rosenbrock_grad_f(x,fx);
	printv(fx);

	// Himmelblau
	printf("Testing if the implementation works on Himmelblau.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton_num(Himmelblau_grad_f,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Himmelblau_grad_f(x,fx);
	printv(fx);

	// part B


	gsl_vector_free(x);
	gsl_vector_free(fx);

	return EXIT_SUCCESS;
}