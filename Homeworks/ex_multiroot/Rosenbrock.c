#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <stdio.h>

int 
rosenbrock_grad (const gsl_vector * x, void *params,
              gsl_vector * f)
{
	(void) params; 

	const double x0 = gsl_vector_get (x, 0);
  	const double x1 = gsl_vector_get (x, 1);

  	const double y0 = -2*(1 - x0)-400*(x1-x0*x0)*x0;
  	const double y1 = 200 * (x1 - x0*x0);

  	gsl_vector_set (f, 0, y0);
  	gsl_vector_set (f, 1, y1);

  	return GSL_SUCCESS;
}

double rosenbrock (double *x, double *y) 
{
	double result = pow((1.0-*x),2)+100.0*pow(*y-pow(*x,2), 2);

	return result;
}
int main() {
	const int size = 2;

	/* all variables are pointers so they can be fed
	to the solvers. */

	/* set the type of solver */
	const gsl_multiroot_fsolver_type *Type;
	Type = gsl_multiroot_fsolver_hybrids;

	/* Make the solver object for the purpose of Type*/
	gsl_multiroot_fsolver *sys;
	sys = gsl_multiroot_fsolver_alloc(Type,size);

	/* initialize the variables x and y presented as a 
	x-vector */
  	gsl_vector *x = gsl_vector_alloc (size);
  	double x_initial[2] = {0.0, 0.0};
  	gsl_vector_set (x, 0, x_initial[0]);
  	gsl_vector_set (x, 1, x_initial[1]);

  	/* defining the multiroot function */
	gsl_multiroot_function f  = {&rosenbrock_grad, size, NULL};
	
	/* finally using the information to initialize fsolver*/
	gsl_multiroot_fsolver_set(sys, &f, x);

	/* the iterator int and status parameter*/
	int i = 0;
	int status;
  printf("Rosenbrock part\n");
	do
    {
      i++;

      /* taking a step*/
      status = gsl_multiroot_fsolver_iterate (sys);

      /* now what are all my variables? */
      double x_val = gsl_vector_get(sys->x, 0);
      double y_val = gsl_vector_get(sys->x, 1);

      double f0_val = gsl_vector_get(sys->f, 0);
      double f1_val = gsl_vector_get(sys->f, 1);

      /* the actual rosenbrock value at the point*/
      double f_rosenbrock =  rosenbrock(&x_val, &y_val);
      
      /* print the value of dependant and independant 
      variables alike*/
      printf("iter = %i \t x = %g \t y = %g \n",i, x_val, y_val);

      /*stop if something was wrong with iteration */
      if (status) break;

      /*check if close enough to zero*/
      status = gsl_multiroot_test_residual (sys->f, 1e-10);
    }
    while (status == GSL_CONTINUE && i < 1000);

  	/*printf ("status = %s\n", gsl_strerror (status));*/

  	gsl_multiroot_fsolver_free (sys);
  	gsl_vector_free (x);

	return 0;
}