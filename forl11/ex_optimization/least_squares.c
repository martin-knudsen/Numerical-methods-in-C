#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>

/*structure of the data inspired by dmitri's implementation*/
struct data {double *t, *y, *e;};

/*the equation used for fitting */
double fit_equation(double t, double A, double B, double T) {
	double f_val = A*exp(-t/T)+B;

	return f_val; 
}

/* the whole least squares method for judging the fit 
which we want to minimize*/
double master(const gsl_vector *x, void *params) {

	/* get data */
	struct data *parameters = (struct data*) params;
	double *t = parameters->t;
	double *y = parameters->y;
	double *e = parameters->e; 

	/* variables to be minimized*/

	double A = gsl_vector_get(x, 0);
	double B = gsl_vector_get(x, 1);
	double T = gsl_vector_get(x, 2);

	int size = sizeof(t);

	double sum = 0; 
	double f_val;

	for(int i=0; i<size; i++) {
		f_val =fit_equation(t[i], A, B, T);

		sum += pow(f_val-y[i],2)/pow(e[i],2);
	}

	return sum;
}

int main() {
	
	/* dimensions of the problem */
	const int size = 3;

	/* experimental data in three arrays*/
	double t_data[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
	double y_data[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
	double e_data[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};

	struct data parameters;
	parameters.t=t_data;
	parameters.y=y_data;
	parameters.e=e_data;


	gsl_multimin_function fun;
	fun.f=&master;
	fun.params=(void*)&parameters;
	fun.n=size;

	gsl_multimin_fminimizer_type *Type \
	= gsl_multimin_fminimizer_nmsimplex2;

	gsl_vector *x = gsl_vector_alloc(3);
	gsl_vector_set(x, 0, 5.0);
	gsl_vector_set(x, 1, 0.0);
	gsl_vector_set(x, 2, 0.5);

	gsl_vector *stepsize = gsl_vector_alloc(3);
	gsl_vector_set_all(stepsize, 1.0);

	gsl_multimin_fminimizer *mini\
	= gsl_multimin_fminimizer_alloc(Type, size);

	gsl_multimin_fminimizer_set(mini, &fun, x, stepsize);

	int status, i=0;
	double it_size;
	double eps_abs = 1e-3;

	double A_val, B_val, T_val, sum_val;

	FILE *my_file;
	my_file = fopen( "least_squares.txt", "w+");

	do {
		i++;

		status = gsl_multimin_fminimizer_iterate(mini);

		if(status) break;

		it_size = gsl_multimin_fminimizer_size(mini);
		status = gsl_multimin_test_size(it_size, eps_abs);

		if(status == GSL_SUCCESS) {
			printf("We have converged to a minimum as seen in the .txt file\n");
		}

		A_val = gsl_vector_get(mini->x, 0);
		B_val = gsl_vector_get(mini->x, 1);
		T_val = gsl_vector_get(mini->x, 2);
		sum_val = mini->fval;

		fprintf(my_file, "%i\t %g\t %g \t %g\t %g\t %g\n"\
			, i, A_val, B_val, T_val, sum_val, it_size);

	}
	while(status ==GSL_CONTINUE && i<500);

	fclose(my_file);

	printf("The lifetime is estimated to be T=%g\n", T_val);

	FILE *plot_file;
	plot_file = fopen("least_squares_plot.txt", "w+");

	int size_data =sizeof(t_data)/sizeof(t_data[0]);

	for(int j=0; j<size_data; j++) {
		double f_value = fit_equation(t_data[j],\
		 A_val, B_val, T_val);

		fprintf(plot_file, "%g \t %g \t %g\n", \
			t_data[j], y_data[j], f_value);
	}
    fclose(plot_file);

    /* for plotting the fitequation*/

    double t = 0.0, tmax = 10.0, tdelta = 0.02;

    FILE *plot_fit_file;
	plot_fit_file = fopen("least_squares_plot_fit.txt", "w+");

    for(t;t<=tmax;t += tdelta) {
    	double fvalue = fit_equation(t, A_val, B_val, T_val);

    	fprintf(plot_fit_file, "%g\t %g \n", t, fvalue);
    }
    fclose(plot_fit_file);

	gsl_vector_free(x);
  	gsl_vector_free(stepsize);
  	gsl_multimin_fminimizer_free (mini);
	return EXIT_SUCCESS;
}