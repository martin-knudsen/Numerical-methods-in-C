#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "ANN.h"
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
	
	printf("Solving y'=y*(1-y), x∈[-5,5], y(0)=0.5 \n");
	
	// number of points for training	
	const int N = 20;  
	
	// interval limits
	double a=-5.0, b=5.0; 

	// make vector of 20 points in interval
	gsl_vector* xvec=gsl_vector_alloc(N);
	double xk;
	for(int k=1; k<=N; k++){
		xk=a+(b-a)*(k-1)/(N-1);
		gsl_vector_set(xvec,k-1,xk);
	}
	//printv(xvec);
	
	printf("Plot can be seen in plot.svg\n");
	// ANN stuff
	// number of hidden neurons
	int number_of_hidden_neurons_log = 7;
	int number_of_hidden_neurons_gauss = 7;
	double activation_function(double x) {
		return x*exp(-x*x);
	}
	double activation_function_der(double x) {
		return exp(-x*x)*(1-2*x*x);
	}

	// start logistic
	double x0_log=0.0, y0_log=0.5, x0_gauss=0.0, y0_gauss=1.0;
	
	double ym_log(double x, double y){ 
		return y*(1-y);
	}

	double y_log_analytic(double x){ 
		double C = exp(x0_log)/y0_log-exp(x0_log);
		return exp(x)/(exp(x)+C);
	}
	double ym_gauss(double x, double y) {
		return -x*y;
	}

	double y_gauss_analytic(double x){ 
		double C = y0_gauss/exp(-0.5*x0_gauss*x0_gauss);
		return C*exp(-0.5*x*x);
	}

	
	ann* network1 = ann_alloc(number_of_hidden_neurons_log, &activation_function, &activation_function_der);
	for(int i=0; i<number_of_hidden_neurons_log;i++){
		gsl_vector_set(network1->data,i*3,a+(b-a)/(number_of_hidden_neurons_log-1)*i);
		gsl_vector_set(network1->data,i*3+1,1);
		gsl_vector_set(network1->data,i*3+2,1);

	}
	ann* network2 = ann_alloc(number_of_hidden_neurons_gauss, &activation_function, &activation_function_der);
	for(int i=0; i<number_of_hidden_neurons_gauss;i++){
		gsl_vector_set(network2->data,i*3,a+(b-a)/number_of_hidden_neurons_gauss*i);
		gsl_vector_set(network2->data,i*3+1,1);
		gsl_vector_set(network2->data,i*3+2,1);

	}
	//printv(network->data);
	ann_train(network1,xvec,&ym_log,x0_log,y0_log);
	ann_train(network2,xvec,&ym_gauss,x0_gauss,y0_gauss);
	// ANN stuff end
	FILE *logdata = fopen("logdata.txt", "w+");
	FILE *gaussdata = fopen("gaussdata.txt", "w+");

	double delta_x = 10.0/100.0, x_start = -5.0, x_slut=5.0, x_it;
	double ANN_result1, analytic_result1,ANN_result2, analytic_result2, deri;
	for(x_it = x_start; x_it <x_slut; x_it += delta_x) {
		ann_feed_forward(network1,x_it,&ANN_result1,&deri);
		ann_feed_forward(network2,x_it,&ANN_result2,&deri);
		analytic_result1 = y_log_analytic(x_it);
		analytic_result2 = y_gauss_analytic(x_it);
		fprintf(logdata, "%g \t %g \t %g\n", x_it, ANN_result1, analytic_result1);
		fprintf(gaussdata, "%g \t %g \t %g\n", x_it, ANN_result2, analytic_result2);
	}

	fclose(logdata);
	fclose(gaussdata);

	gsl_vector_free(xvec);
	ann_free(network1);
	ann_free(network2);
	return EXIT_SUCCESS;
}