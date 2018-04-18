#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
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
	int function_calls = 0;
	bool ANALYTIC = false;


	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_matrix* J = gsl_matrix_alloc(n,n);

	void system_f(gsl_vector* v, gsl_vector* fx, gsl_matrix* J){
	double A = 10000.0;
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx0 = A*x*y -1; 
	double fx1 = exp(-x) + exp(-y) -(1.0 + 1.0/A);
	double J00 = A*y;
	double J01 = A*x;
	double J10 = -exp(-x);
	double J11 = -exp(-y);
	gsl_matrix_set(J,0,0,J00);
	gsl_matrix_set(J,0,1,J01);
	gsl_matrix_set(J,1,0,J10);
	gsl_matrix_set(J,1,1,J11);
	gsl_vector_set(fx,0,fx0);
	gsl_vector_set(fx,1,fx1);
	function_calls++;
	}

	void Rosenbrock_grad_f(gsl_vector* v, gsl_vector* fx, gsl_matrix* J){
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx0 = 2*x-2+400*(x*x*x-y*x); 
	double fx1 = 200*(y-x*x);
	double J00 = 2+100*(12*x*x-4*y);
	double J01 = -400*x;
	double J10 = -400*x;
	double J11 = 200;
	gsl_matrix_set(J,0,0,J00);
	gsl_matrix_set(J,0,1,J01);
	gsl_matrix_set(J,1,0,J10);
	gsl_matrix_set(J,1,1,J11);
	gsl_vector_set(fx,0,fx0);
	gsl_vector_set(fx,1,fx1);
	function_calls++;
	}

	void Himmelblau_grad_f(gsl_vector* v, gsl_vector* fx, gsl_matrix* J){
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx0 = 4*x*(x*x+y-11)+2*(x+y*y-7); 
	double fx1 = 2*(x*x+y-11)+4*y*(x+y*y-7);
	double J00 = 4*(3*x*x+y-11)+2;
	double J01 = 4*(x+y);
	double J10 = 4*(x+y);
	double J11 = 2+4*(x+3*y*y-7);
	gsl_matrix_set(J,0,0,J00);
	gsl_matrix_set(J,0,1,J01);
	gsl_matrix_set(J,1,0,J10);
	gsl_matrix_set(J,1,1,J11);
	gsl_vector_set(fx,0,fx0);
	gsl_vector_set(fx,1,fx1);
	function_calls++;
	}
	
	
	// system
	printf("A. Numerical Jacobian\n");
	double dx=1e-9, epsilon=1e-3;
	printf("Testing if the implementation works on the system of equations provided.\n");
	
	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton(system_f,x,dx,epsilon, ANALYTIC);	
	
	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is (should be close to zero-vector)\n");
	system_f(x,fx,J);
	printv(fx);
	printf("Number of function calls: %i\n",function_calls);
	printf("\n");
	// Rosenbrock
	function_calls = 0;
	dx=1e-6, epsilon=1e-3;
	printf("Testing if the implementation works on Rosenbrock.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton(Rosenbrock_grad_f,x,dx,epsilon, ANALYTIC);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Rosenbrock_grad_f(x,fx,J);
	printv(fx);
	printf("Number of function calls: %i\n",function_calls);
	printf("\n");

	// Himmelblau
	function_calls = 0;
	printf("Testing if the implementation works on Himmelblau.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton(Himmelblau_grad_f,x,dx,epsilon, ANALYTIC);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Himmelblau_grad_f(x,fx,J);
	printv(fx);
	printf("Number of function calls: %i\n",function_calls);

	// part B
	printf("\n\n");
	printf("B. using analytical Jacobian\n");
	ANALYTIC = true;

	// system
	function_calls = 0;
	dx=1e-6, epsilon=1e-3;
	printf("Testing if the implementation works on the system of equations provided.\n");
	
	gsl_vector_set(x,0,1);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton(system_f,x,dx,epsilon, ANALYTIC);	
	
	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is (should be close to zero-vector)\n");
	system_f(x,fx,J);
	printv(fx);
	printf("Number of function calls: %i\n",function_calls);
	printf("\n");

	// Rosenbrock
	function_calls = 0;
	dx=1e-8, epsilon=1e-3;
	printf("Testing if the implementation works on Rosenbrock.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton(Rosenbrock_grad_f,x,dx,epsilon, ANALYTIC);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Rosenbrock_grad_f(x,fx,J);
	printv(fx);
	printf("Number of function calls: %i\n",function_calls);
	printf("\n");

	// Himmelblau
	function_calls = 0;
	printf("Testing if the implementation works on Himmelblau.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton(Himmelblau_grad_f,x,dx,epsilon, ANALYTIC);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Himmelblau_grad_f(x,fx,J);
	printv(fx);
	printf("Number of function calls: %i\n",function_calls);
	
	printf("\n");
	printf("As one can see the use of an analytical jacobian\n
		does not have a dramatic effect on the stepsize\n
		and even increases it a bit in the Rosenbrock function. \n
		This could be a coincidence as a slightly wrong jacobian\n
		might luckily force the jacobian in a good direction of\n
		function space maybe to not get stuck in a bad local place\n
		Probably if we would have more variables it would matter more\n
		However the analytical version is superior everywhere in terms\n
		of speed as it limits the number of calls to the function significantly\n
		which it uses 2-3 times less.
		");


	printf("\n\n");


	//part C
	printf("C. Advanced backtracking");

	gsl_vector_free(x);
	gsl_vector_free(fx);
	gsl_matrix_free(J);

	return EXIT_SUCCESS;
}