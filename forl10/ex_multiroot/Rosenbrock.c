#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <stdio.h>

int 
rosenbrock (const gsl_vector * x, void *params,
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

int main() {
	int size = 2;

	const gsl_multiroot_fsolver_type *Type;
	Type = gsl_multiroot_fsolver_hybrids;

	gsl_multiroot_fsolver *sys;
	sys = gsl_multiroot_fsolver_alloc(Type,size);

  	gsl_vector *x = gsl_vector_alloc (size);
  	double x_initial[2] = {0.0, 0.0};

  	gsl_vector_set (x, 0, x_initial[0]);
  	gsl_vector_set (x, 1, x_initial[1]);

	gsl_multiroot_function f  = {&rosenbrock, size, NULL};
	
	gsl_multiroot_fsolver_set(sys, &f, x);

	int i = 0;
	int status;

	do
    {
      i++;
      status = gsl_multiroot_fsolver_iterate (sys);

      double x_val = gsl_vector_get(sys->x, 0);
      double y_val = gsl_vector_get(sys->x, 1);

      double f0_val = gsl_vector_get(sys->f, 0);
      double f1_val = gsl_vector_get(sys->f, 1);

      printf("i = %i \t x=%g \t y=%g \t f0 = %g \t \
       f1 = %g\n",i, x_val, y_val, f0_val, f1_val);

      if (status)   
        break;

      status = gsl_multiroot_test_residual (sys->f, 1e-10);
    }
  while (status == GSL_CONTINUE && i < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (sys);
  gsl_vector_free (x);

	return 0;
}