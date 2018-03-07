#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h> 
#include <stdlib.h>

int 
diff_eq (double r, const double y[], double dydt[], void * params) {
  
  double E = *(double *)params;
  dydt[0] = y[1];
  dydt[1] = -2*(E*y[0]+1/r*y[0]);

  return GSL_SUCCESS;
}


 
double F_E (double E, double r)            
{
  const double rmin = 1e-4;
  if(r<rmin) return r-r*r;

	gsl_odeiv2_system sys;
  sys.function = diff_eq;
  sys.jacobian = NULL;
  sys.dimension = 2;
  sys.params = (void*)&E;

  double hstart = 1e-3; double epsabs = 1e-6; double epsrel = 1e-6; 

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);

  double t = rmin; 
  double y[2] = {t-t*t, 1-2*t};
  /*
  printf("-----------------------------------------------\n");
  printf("Now solving the equatorial equation for eps = 0\n");
  printf("Resulting in the differential equation u''=1-u\n");
  */
      
  int status = gsl_odeiv2_driver_apply (d, &t, r, y);
  if (status != GSL_SUCCESS)
  {
    printf ("error, return value=%i\n", status);
    
  }

  gsl_odeiv2_driver_free (d);

  return y[0];
}

int 
F_E_control (const gsl_vector *x, void *params, gsl_vector *f)
{
  double E = gsl_vector_get(x, 0);

  double rmax =*(double*)params;

  double f_val =F_E(E, rmax);

  gsl_vector_set(f, 0, f_val);

  return GSL_SUCCESS;
}

double exact_F(double r) {
  double result = r*exp(-r);

  return result;
}

int main() {
	const int size = 1;
  const double rmax = 8.0; 
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
  	
    double x_initial = -13.6;
  	gsl_vector_set (x, 0, x_initial);
  
  	/* defining the multiroot function */
	gsl_multiroot_function f  = {&F_E_control, size, (void*)&rmax};
	
	/* finally using the information to initialize fsolver*/
	gsl_multiroot_fsolver_set(sys, &f, x);

	/* the iterator int and status parameter*/
	int i = 0;
	int status;
  double E_val = 0.0, f0_val = 0.0; 
	do
    {
      i++;

      /* taking a step*/
      status = gsl_multiroot_fsolver_iterate (sys);

      /* now what are all my variables? */
      E_val = gsl_vector_get(sys->x, 0);
      f0_val = gsl_vector_get(sys->f, 0);
      
      /* print the value of dependant and independant 
      variables alike*/
      printf("iter = %i \t E = %g \t  M(E) = %g \n",i, E_val, f0_val);
      /*stop if something was wrong with iteration */
      if (status) 
      {
      printf ("status = %s\n", gsl_strerror (status));
      break;
      }

      /*check if close enough to zero*/
      status = gsl_multiroot_test_residual (sys->f, 1e-10);
    }
    while (status == GSL_CONTINUE && i < 1000);

  	printf ("status = %s\n", gsl_strerror (status));
    printf("So the solution for rmax = %g is E=%g\n", rmax, E_val);
  	
    FILE *my_file;
    my_file = fopen("shooting.txt","w+");

    double rmin = 1e-4;
    for(double r=rmin; r<=rmax; r+=0.02){
      fprintf(my_file, "%g \t %g \t %g\n", r, F_E(E_val, r), exact_F(r));
    }

    fclose(my_file);


    gsl_multiroot_fsolver_free (sys);
  	gsl_vector_free (x);

	return 0;
}