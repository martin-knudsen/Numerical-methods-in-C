#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Newton.h"
#include "downhill_simplex.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.8f" //format of print "7 width, 3 digits after comma" 


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
	int n=2;
	int function_calls = 0;
	double epsilon=1e-3, dx=1e-6;

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

	double master_fit(gsl_vector* v){
	double A = gsl_vector_get(v,0);
	double B =	gsl_vector_get(v,1);
	double T =	gsl_vector_get(v,2);
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);
	double fx = 0;
	for(int i=0;i<N;i++){
		fx +=pow((A*exp(-t[i]/T)+B-y[i])/e[i],2);
	}
	function_calls++;
	return fx;
	}

	double Rosenbrock_f_downhill(double* v){
	double x = v[0];
	double y = v[1];
	double fx = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
	function_calls++;
	return fx;
	}

	double Himmelblau_f_downhill(double* v){
	double x = v[0];
	double y = v[1];
	double fx = (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
	function_calls++;
	return fx;
	}

	double master_fit_downhill(double* v){
	double A = v[0];
	double B = v[1];
	double T = v[2];
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);
	double fx = 0;
	for(int i=0;i<N;i++){
		fx +=pow((A*exp(-t[i]/T)+B-y[i])/e[i],2);
	}
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
	printf("Here df(x) is\n");
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
	printf("Here df(x) is\n");
	Himmelblau_f(x,df,H);
	printv(df);
	printf("Number of function calls: %i\n",function_calls);
	printf("We can see that this algorithm is much less efficient as\n"
			"the root finding one. Again this might be because of parameters\n"
			"the bigger dx. I have set alfa=0.51 it seems good for both functions \n");
	printf("\n\n");

	printf("B. Quasi Newton with Broyden analytical gradient\n");
	printf("i)\n");
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
	printf("Here df(x) is\n");
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

	newton_minimization_Broyden(Himmelblau_f_Broyden,x,dx,epsilon);	

	printf("The root found is: \n");
	printv(x);
	printf("Here df(x) is\n");
	Himmelblau_f(x,df,H);
	printv(df);
	printf("Number of function calls: %i\n",function_calls);
	printf("Something very interesting happened when using Broyden:\n"
			"The minima for Himmelblau found is a different one, but still correct\n"
			"since Himmelblau has several minima.\n");

	printf("\n");
	printf("iii) Comparison\n"
			"If we compare explicit hessian to Broyden, we see that \n"
			"Explicit Hessian uses many fewer steps than updating without Hessian\n"
			"This does no nessesarily mean that classical Newton is better\n"
			"actually possibly broyden is still faster because it doesn't \n"
			"solve the linear system of equations to find delta x at each step and doesn't\n"
			"calculate the hessian each time the function is run. \n"
			"If we compare to the root finding algorithms they take even less steps\n"
			"and seem to be the most effective overall.\n"
			"Again however it is difficult to say overall because it might be a parameter issue\n");

	n=3;
	printf("\n");
	printf("iv) Testing at fit problem\n");
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);

	gsl_vector* x_fit = gsl_vector_alloc(n);
	gsl_vector* df_fit = gsl_vector_alloc(n);

	FILE* data = fopen("data.txt","w+");
	for(int i=0; i<N; i++){
		fprintf(data,"%g\t%g\t%g\n",t[i],y[i],e[i]);
	}
	fclose(data);

	printf("I have found the gradient of the master funciton numerically\n");
	function_calls = 0;
	gsl_vector_set(x_fit,0,5);
	gsl_vector_set(x_fit,1,0);
	gsl_vector_set(x_fit,2,5);
	
	printf("I set dx=%g and epsilon=%g\n",dx, epsilon);
	printf("xstart is (qualified guess from data):\n");
	printv(x_fit);

	newton_minimization_Broyden_num_gradient(master_fit,x_fit,dx,epsilon);	

	printf("The root found is: \n");
	printv(x_fit);
	printf("Here df(x) is\n");
	gradient_num(master_fit,x_fit,df_fit,dx);
	printv(df_fit);
	printf("Number of function calls: %i\n",function_calls);

	double tstart=0.0, tslut=10.0, deltat=tslut/200;
	double A=gsl_vector_get(x_fit,0);
	double B=gsl_vector_get(x_fit,1);
	double T=gsl_vector_get(x_fit,2);
	double ft;
	FILE* fitdata = fopen("fit.txt","w+");

	for(double t=0;t<tslut;t+=deltat){
		ft=A*exp(-t/T)+B;
		fprintf(fitdata, "%g\t%g\n",t, ft);
	}
	fclose(fitdata);
	printf("\n\n");
	


	// C
	
	printf("C. downhill simplex method\n");
	n=2;
	double** simplex= (double **) malloc ((n+1) * sizeof(double *));
	for (int i=0;i<=n;i++) {
    simplex[i] = (double *) malloc (n * sizeof(double));
	}
	// Rosenbrock
	
	function_calls = 0;
	printf("Testing if the implementation works on Rosenbrock.\n");
	simplex[0][0]=0.0;
	simplex[0][1]=1.0;
	simplex[1][0]=0.34;
	simplex[1][1]=5.0;
	simplex[2][0]=-10;
	simplex[2][1]=0.0;
	printf("I set epsilon=%g\n",epsilon);
	printf("I set the simplex to some \"random\" values:\n");
	for(int i=0; i<n+1;i++){
		for(int j=0; j<n;j++){printf("%g\t",simplex[i][j]);}
		printf("\n");
	}
	int k= downhill_simplexx(Rosenbrock_f_downhill,simplex,n,epsilon);
	
	printf("The root found is: \n");
	for(int i=0; i<n;i++){
		printf("%g\t",simplex[0][i]);
	}
	printf("\n");
	printf("The number of steps: %i\n",k);
	printf("Number of function calls: %i\n",function_calls);
	printf("\n");
	
	// Himmelblau
	
	function_calls = 0;
	printf("Testing if the implementation works on Himmelblau.\n");
	simplex[0][0]=0.0;
	simplex[0][1]=1.0;
	simplex[1][0]=0.34;
	simplex[1][1]=5.0;
	simplex[2][0]=-10;
	simplex[2][1]=0.0;
	printf("I set epsilon=%g\n",epsilon);
	printf("I set the simplex to some \"random\" values:\n");
	for(int i=0; i<n+1;i++){
		for(int j=0; j<n;j++){printf("%g\t",simplex[i][j]);}
		printf("\n");
	}
	k= downhill_simplexx(Himmelblau_f_downhill,simplex,n,epsilon);
	
	printf("The root found is: \n");
	for(int i=0; i<n;i++){
		printf("%g\t",simplex[0][i]);
	}
	printf("\n");
	printf("The number of steps: %i\n",k);
	printf("Number of function calls: %i\n",function_calls);
	printf("Again we have another value for the minimum, but it's \n"
		"just another of the 5 minimums of the function!\n");
	printf("\n");

	// Master fit
	n=3;
	function_calls = 0;
	printf("Testing if the implementation works on the fit.\n");
	double** simplex_fit= (double **) malloc ((n+1) * sizeof(double *));
	for (int i=0;i<=n;i++) {
    simplex_fit[i] = (double *) malloc (n * sizeof(double));
	}

	simplex_fit[0][0]=5.0;
	simplex_fit[0][1]=0.0;
	simplex_fit[0][2]=5.0;
	simplex_fit[1][0]=10.0;
	simplex_fit[1][1]=0.1;
	simplex_fit[1][2]=10.0;
	simplex_fit[2][0]=3;
	simplex_fit[2][1]=0.5;
	simplex_fit[2][2]=5.0;
	simplex_fit[3][0]=7;
	simplex_fit[3][1]=0.4;
	simplex_fit[3][2]=3;
	printf("I set epsilon=%g\n",epsilon);
	printf("I set the simplex to some \"random\" values:\n");
	for(int i=0; i<n+1;i++){
		for(int j=0; j<n;j++){printf("%g\t",simplex_fit[i][j]);}
		printf("\n");
	}
	k= downhill_simplexx(master_fit_downhill,simplex_fit,n,epsilon);
	
	printf("The root found is: \n");
	for(int i=0; i<n;i++){
		printf("%g\t",simplex_fit[0][i]);
	}
	printf("\n");
	printf("The number of steps: %i\n",k);
	printf("Number of function calls: %i\n",function_calls);
	gsl_vector_set(x_fit,0,simplex_fit[0][0]);
	gsl_vector_set(x_fit,1,simplex_fit[0][1]);
	gsl_vector_set(x_fit,2,simplex_fit[0][2]);

	gradient_num(master_fit,x_fit,df_fit,dx);
	printf("Here df(x) is\n");
	printv(df_fit);
	printf("So it is almost the same solution as before \n");
	printf("\n");

	gsl_vector_free(x);
	gsl_vector_free(df);
	gsl_vector_free(x_fit);
	gsl_vector_free(df_fit);
	gsl_matrix_free(H);		

	for (int i=0;i<=2;i++) {
    	free (simplex[i]);
  	}
	free(simplex);

	for (int i=0;i<=3;i++) {
    	free (simplex_fit[i]);
  	}
	free(simplex_fit);


	return EXIT_SUCCESS;
}