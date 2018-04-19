#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Newton.h"
#include <gsl/gsl_multiroots.h>
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
	double epsilon=1e-3, dx=1e-3;

	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_matrix* H = gsl_matrix_alloc(n,n);

	double Rosenbrock_f(gsl_vector* v, gsl_vector* df, gsl_matrix* H){
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
	double df0 = 2*x-2+400*(x*x*x-y*x); 
	double df1 = 200*(y-x*x);
	double H00 = 2+100*(12*x*x-4*y);
	double H01 = -400*x;
	double H10 = -400*x;
	double H11 = 200;
	gsl_matrix_set(H,0,0,H00);
	gsl_matrix_set(H,0,1,H01);
	gsl_matrix_set(H,1,0,H10);
	gsl_matrix_set(H,1,1,H11);
	gsl_vector_set(df,0,df0);
	gsl_vector_set(df,1,df1);
	function_calls++;
	return fx;
	}

	double Himmelblau_f(gsl_vector* v, gsl_vector* df, gsl_matrix* H){
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx = (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
	double df0 = 4*x*(x*x+y-11)+2*(x+y*y-7); 
	double df1 = 2*(x*x+y-11)+4*y*(x+y*y-7);
	double H00 = 4*(3*x*x+y-11)+2;
	double H01 = 4*(x+y);
	double H10 = 4*(x+y);
	double H11 = 2+4*(x+3*y*y-7);
	gsl_matrix_set(H,0,0,H00);
	gsl_matrix_set(H,0,1,H01);
	gsl_matrix_set(H,1,0,H10);
	gsl_matrix_set(H,1,1,H11);
	gsl_vector_set(df,0,df0);
	gsl_vector_set(df,1,df1);
	function_calls++;
	return fx;
	}

	double Rosenbrock_f_Broyden(gsl_vector* v, gsl_vector* df){
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
	double df0 = 2*x-2+400*(x*x*x-y*x); 
	double df1 = 200*(y-x*x);
	gsl_vector_set(df,0,df0);
	gsl_vector_set(df,1,df1);
	function_calls++;
	return fx;
	}

	double Himmelblau_f_Broyden(gsl_vector* v, gsl_vector* df){
	double x = gsl_vector_get(v,0);
	double y =	gsl_vector_get(v,1);
	double fx = (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
	double df0 = 4*x*(x*x+y-11)+2*(x+y*y-7); 
	double df1 = 2*(x*x+y-11)+4*y*(x+y*y-7);
	gsl_vector_set(df,0,df0);
	gsl_vector_set(df,1,df1);
	function_calls++;
	return fx;
	}
	
	printf("A. Newton minimization\n");
	
	// Rosenbrock
	function_calls = 0;
	printf("Testing if the implementation works on Rosenbrock.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton_minimization(Rosenbrock_f,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Rosenbrock_f(x,df,H);
	printv(df);
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

	newton_minimization(Himmelblau_f,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Himmelblau_f(x,df,H);
	printv(df);
	printf("Number of function calls: %i\n",function_calls);
	printf("We can see that this algorithm is much less efficient as\n"
			"the root finding one. Again this might be because of parameters\n"
			"the bigger dx. I have set alfa=0.51 it seems good for both functions \n");
	printf("\n\n");

	printf("B. Quasi Newton with Broyden analytical gradient\n");
	
	// Rosenbrock
	function_calls = 0;
	printf("Testing if the implementation works on Rosenbrock.\n");
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,0);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is:\n");
	printv(x);

	newton_minimization_Broyden(Rosenbrock_f_Broyden,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Rosenbrock_f(x,df,H);
	printv(df);
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

	newton_minimizationBroyden(Himmelblau_f_Broyden,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here f(x) is\n");
	Himmelblau_f_B(x,df,H);
	printv(df);
	printf("Number of function calls: %i\n",function_calls);
	gsl_vector_free(x);
	gsl_vector_free(df);
	gsl_matrix_free(H);		

	return EXIT_SUCCESS;
}