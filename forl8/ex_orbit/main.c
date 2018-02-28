#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

int 
function (double t, const double y[], double dydt[], void * params) {
	
	(void)(t); 
	(void)(params); 
	dydt[0] = y[0]*(1-y[0]);

	return GSL_SUCCESS;
}
int 
jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params) {
	
	(void)(t); 
	(void)(params); 
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);
  	gsl_matrix * m = &dfdy_mat.matrix;
  	gsl_matrix_set(m, 1, 1, 1-2.0*y[0]);
  	dfdt[0]=0.0;
	return GSL_SUCCESS;
}

int main() {


	gsl_odeiv2_system sys = {function, NULL, 1, NULL};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

	double x = 0.0; double x1 = 3.0;
	double y[1]= {0.5};

	printf("Now testing the solution for y'=y(1-y)\n");
	printf("by using the gsl_odeiv2_step_rkf45 stepper");
	
	for (int i = 1; i<=100 ; i++) {
		double xi = i * x1 /100.0;
		int status = gsl_odeiv2_driver_apply (d, &x, xi, y);
		if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      printf ("%.5e %.5e \n", x, y[0]);

	}
	gsl_odeiv2_driver_free (d);

	double anal_solution = exp(3)/(exp(3)+1);

	printf("The analytical solution at x=3 is %g\n", anal_solution);
	return 0;
}