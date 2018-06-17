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
	gsl_vector* x=gsl_vector_alloc(N);
	for(int k=1; k<=N; k++){
		xk=a+(b-a)*(k-1)/(N-1)
		gsl_vector_set(x,k-1,xk);
	}
	//printv(xvec);
	
	printf("Plot can be seen in plot.svg\n");
	
	// ANN stuff
	// number of hidden neurons
	int number_of_hidden_neurons = 5;
	double activation_function(double x) {
		return x*exp(-x*x);
	}

	double activation_function_der(double x) {
		return exp(-x*x)*(1-2*x*x);
	}

	// start logistic

	double ym_log(double x, double y) return y*(1-y);
	double ym_gauss(double x, double y) return -x*y;

	double x0_log=0.0, y0_log=0.5, x0_gauss=0.0, y0_gauss=1.0;

	ann* network1 = ann_alloc(number_of_hidden_neurons, &activation_function, &activation_function_der);
	for(int i=0; i<number_of_hidden_neurons;i++){
		gsl_vector_set(network->data,i*3,2*M_PI*i/number_of_hidden_neurons);
		gsl_vector_set(network->data,i*3+1,1);
		gsl_vector_set(network->data,i*3+2,1);

	}

	ann* network2 = ann_alloc(number_of_hidden_neurons, &activation_function, &activation_function_der);
	for(int i=0; i<number_of_hidden_neurons;i++){
		gsl_vector_set(network->data,i*3,2*M_PI*i/number_of_hidden_neurons);
		gsl_vector_set(network->data,i*3+1,1);
		gsl_vector_set(network->data,i*3+2,1);

	}
	//printv(network->data);

	ann_train(network1,xlist,&ym_log,x0_log,y0_log);
	//ann_train(network2,xlist,&ym_gauss,x0_gauss,y0_gauss);
	// ANN stuff end

	FILE *lineardata = fopen("lineardata.txt", "w+");

	/*
	double delta_x = 2*M_PI/100, x_start = -5.0, x_slut=5.0, x_it;
	double spline_result;
	for(x_it = x_start; x_it <x_slut; x_it += delta_x) {
		spline_result = ann_feed_forward(network,x_it);
		//printf("%g\n",ann_feed_forward(network,x_it));
		fprintf(lineardata, "%g \t %g \n", x_it, spline_result);
	}

	fclose(lineardata);

	FILE *data_points  =fopen("data_points.txt", "w+");
	for(int i=0; i<n; i++) {
		fprintf(data_points, "%g \t %g\n", x[i], y[i]);
	}
	fclose(data_points);
	*/

	gsl_vector_free(xvec);
	ann_free(network1);
	ann_free(network2);
	return EXIT_SUCCESS;
}