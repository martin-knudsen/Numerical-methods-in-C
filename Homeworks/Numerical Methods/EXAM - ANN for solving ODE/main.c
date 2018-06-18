#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "ANN.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.6f" //format of print "7 width, 3 digits after comma" 

void printv(gsl_vector *A){
	for(int i=0;i<A->size;i++){
		printf(FMT,gsl_vector_get(A,i));
		printf("\n");
	}
}

int main(void) {
	
	printf("The exam problem solved is 'Artificial neural network (ANN) for solving ODE'. \n");
	printf("Solving logistic function\n"
		"y'=y*(1-y), x∈[-5,5], y(0)=0.5. \n");
	printf("Solving gaussian function\n"
		"y'=-x*y, x∈[-5,5], y(0)=1. \n");

	
	// number of points for training	
	const int N = 25;  
	
	// interval limits
	double a=-5.0, b=5.0; 

	// Initial conditions 
	double x0_log=0.0, y0_log=0.5, x0_gauss=0.0, y0_gauss=1.0;
	
	// make vector of N points in interval
	gsl_vector* xvec=gsl_vector_alloc(N);
	double xk;
	for(int k=1; k<=N; k++){
		xk=a+(b-a)*(k-1)/(N-1);
		gsl_vector_set(xvec,k-1,xk);
	}
	
	int number_of_hidden_neurons_log = 6;
	int number_of_hidden_neurons_gauss = 6;
	printf("Number of hidden neurons for both=%i.\n",number_of_hidden_neurons_log);
	printf("Number of points in the interval used for training=%i\n",N);

	printf("Activation function is gaussian wavelet for both.\n");
	printf("The analytical derivative of the activation function is used.\n"); 
	double activation_function(double x) {
		return x*exp(-x*x);
	}
	double activation_function_der(double x) {
		return exp(-x*x)*(1-2*x*x);
	}

	// logarithmic function
	double ym_log(double x, double y){ 
		return y*(1-y);
	}

	// its analytical solution found by seperation of variables
	double y_log_analytic(double x){ 
		double C = exp(x0_log)/y0_log-exp(x0_log);
		return exp(x)/(exp(x)+C);
	}

	// gaussian function
	double ym_gauss(double x, double y) {
		return -x*y;
	}

	// analytical solution found by seperation of variables again
	double y_gauss_analytic(double x){ 
		double C = y0_gauss/exp(-0.5*x0_gauss*x0_gauss);
		return C*exp(-0.5*x*x);
	}

	// allocate the ANN for the logistic function. Now it contains the derivative of the activation function
	// starting with weights of 1 and making sure the offset ai is equally distibuted in the interval
	ann* network_log = ann_alloc(number_of_hidden_neurons_log, &activation_function, &activation_function_der);
	for(int i=0; i<number_of_hidden_neurons_log;i++){
		gsl_vector_set(network_log->data,i*3,a+(b-a)/(number_of_hidden_neurons_log-1)*i);
		gsl_vector_set(network_log->data,i*3+1,1);
		gsl_vector_set(network_log->data,i*3+2,1);
	}

	// do the same for gauss
	ann* network_gauss = ann_alloc(number_of_hidden_neurons_gauss, &activation_function, &activation_function_der);
	for(int i=0; i<number_of_hidden_neurons_gauss;i++){
		gsl_vector_set(network_gauss->data,i*3,a+(b-a)/(number_of_hidden_neurons_gauss-1)*i);
		gsl_vector_set(network_gauss->data,i*3+1,1);
		gsl_vector_set(network_gauss->data,i*3+2,1);
	}

	printf("Now both networks are trained by using my quasi-Newton method\n"
	 "with broyden update implementation from the minimization homework.\n");
	ann_train(network_log,xvec,&ym_log,x0_log,y0_log);
	ann_train(network_gauss,xvec,&ym_gauss,x0_gauss,y0_gauss);

	// this is for printing out data for the plots. 
	

	// Compare 100 points in the interval between the analytical and the ANN solution for both
	// functions. 
	FILE *logdata = fopen("logdata.txt", "w+");
	FILE *gaussdata = fopen("gaussdata.txt", "w+");
	double delta_x = 10.0/100.0, x_start = -5.0, x_slut=5.0, x_it;
	double ANN_result1, analytic_result1,ANN_result2, analytic_result2, deri;
	for(x_it = x_start; x_it <x_slut; x_it += delta_x) {
		ann_feed_forward(network_log,x_it,&ANN_result1,&deri);
		ann_feed_forward(network_gauss,x_it,&ANN_result2,&deri);
		analytic_result1 = y_log_analytic(x_it);
		analytic_result2 = y_gauss_analytic(x_it);
		fprintf(logdata, "%g \t %g \t %g\n", x_it, ANN_result1, analytic_result1);
		fprintf(gaussdata, "%g \t %g \t %g\n", x_it, ANN_result2, analytic_result2);
	}
	fclose(logdata);
	fclose(gaussdata);

	printf("For comparison of how well the ANN solves the equations\n"
		"both the analytical solution found by seperation of variables and the ANN solution\n"
		"for 100 points in the interval are found.\n"
		"This is plotted in plot_log.svg and plot_gauss.svg respectively.\n"
		"As on can see the solutions are approximated very well with exception of\n"
		"the ends of the intervals.\n ");

	gsl_vector_free(xvec);
	ann_free(network_log);
	ann_free(network_gauss);
	return EXIT_SUCCESS;
}