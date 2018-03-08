#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multimin.h>


/* the function to be minimized */
double rosenbrock(const gsl_vector *x, void *params) {
	double *parameters = (double *)params;
	double a = parameters[0];
	double b = parameters[1];

	double x_val = gsl_vector_get(x, 0); 
	double y_val = gsl_vector_get(x, 1);

	double f_val  = pow(a-x_val, 2) + \
	b*pow(y_val-pow(x_val, 2), 2);

	return f_val;
}


int main() {
	
	const int size = 2;
	double parameters[2] ={1.0, 100.0}; 

	/* initializing the multimin function */
	gsl_multimin_function fun;
	fun.f = &rosenbrock;
	fun.n = size;
	fun.params = (void *)parameters;

	/* initializing the type of solver */

	gsl_multimin_fminimizer_type *Type;
	Type = gsl_multimin_fminimizer_nmsimplex2;

	/* init the minimizer itself */
	gsl_multimin_fminimizer *mini;
	mini = gsl_multimin_fminimizer_alloc(Type, size);

	/* set the vectors containing stepsize and variables*/
	gsl_vector *x = gsl_vector_alloc(2);
	gsl_vector_set_all(x, 0);

	gsl_vector *stepsize = gsl_vector_alloc(2);
	gsl_vector_set_all(stepsize, 0.1);

	/* set the minimizer to the function */
	gsl_multimin_fminimizer_set(mini, &fun, x, stepsize);

	/*iterating as usual*/
	int status, i=0; 
	double it_size;
	double eps_abs = 1e-3;

	double x_val, y_val, f_val;

	FILE *my_file;
	my_file = fopen( "rosenbrock_minimum.txt", "w+");
	do {
		i++;
		status = gsl_multimin_fminimizer_iterate(mini);

		if(status) break;

		it_size = gsl_multimin_fminimizer_size(mini);
		status = gsl_multimin_test_size(it_size, eps_abs);

		if(status == GSL_SUCCESS) {
			printf("We have converged to a minimum at \n");
		}

		x_val = gsl_vector_get(mini->x, 0);
		y_val = gsl_vector_get(mini->x, 1);
		f_val = mini->fval;
		


		fprintf(my_file, "%i\t %g\t %g\t %g\t %g\n"\
			, i, x_val, y_val, f_val, it_size);
	}
	while(status == GSL_CONTINUE && i<500);

	fclose(my_file);
	gsl_vector_free(x);
  	gsl_vector_free(stepsize);
  	gsl_multimin_fminimizer_free (mini);

  	return EXIT_SUCCESS;
}