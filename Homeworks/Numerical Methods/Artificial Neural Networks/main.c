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
	
	printf("A. Testing the ANN on the function I used in the interpolant exercise\n");
	
	const int n = 11;  
	double x[11] = {0.0, 2*M_PI/20*2, 2*M_PI/20*4, \
		2*M_PI/20*6, 2*M_PI/20*8, \
		2*M_PI/20*10, 2*M_PI/20*12,  \
		2*M_PI/20*14, 2*M_PI/20*16, \
		2*M_PI/20*18, 2*M_PI/20*20};
	double y[11] = {1.0, 0.80901,  0.30901, \
		 -0.30901,   -0.80901,   -1.0, \
		   -0.80901,    -0.30901,  \
		 0.30901, 0.80901,  1.0}; 

	gsl_vector* yvec=gsl_vector_alloc(n);
	gsl_vector* xvec=gsl_vector_alloc(n);
	for(int i=0; i<n; i++){
		gsl_vector_set(xvec,i,x[i]);
		gsl_vector_set(yvec,i,y[i]);
	}
	//printv(xvec);
	//printv(yvec);
	
	printf("Testing if linear spline works. \n");
	printf("Testfunction is cos(x) from 0 to 2pi. \n");
	printf("Plot can be seen in plot.svg \n");
	
	// ANN stuff
	int number_of_hidden_neurons = 5;
	double activation_function(double x) {
		return x*exp(-x*x);
	}

	ann* network = ann_alloc(number_of_hidden_neurons, &activation_function);
	for(int i=0; i<number_of_hidden_neurons;i++){
		gsl_vector_set(network->data,i*3,2*M_PI*i/number_of_hidden_neurons);
		gsl_vector_set(network->data,i*3+1,1);
		gsl_vector_set(network->data,i*3+2,1);


	}
	//printv(network->data);

	ann_train(network,xvec,yvec);
	// ANN stuff end

	FILE *lineardata = fopen("lineardata.txt", "w+");


	double delta_x = 2*M_PI/100, x_start = 0.0, x_slut=2*M_PI, x_it;
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

	gsl_vector_free(xvec);
	gsl_vector_free(yvec);

	ann_free(network);
	return EXIT_SUCCESS;
}