#ifndef HAVE_Runge_Kutta
#define HAVE_Runge_Kutta

void rkstep12(void f(int n,double x,gsl_vector *y,gsl_vector *dydx),
int n, double x,gsl_vector* yx,double h,gsl_vector* yh,gsl_vector* dy);

int ode_driver(void f(int n,double x,gsl_vector *y,gsl_vector *dydx),
int n,gsl_vector *xlist,gsl_matrix *ylist,
double b,double h,double acc,double eps,int max);

double integrator(double f(double x),
int n,gsl_vector *xlist,gsl_matrix *ylist,
double b,double h,double acc,double eps,int max);

#endif
