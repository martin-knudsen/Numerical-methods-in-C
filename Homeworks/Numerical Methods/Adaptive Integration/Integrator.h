#ifndef HAVE_Integrator
#define HAVE_Integrator

double adapt24(double f(double),double old_f(double),double a,double b,
	double acc,double eps,double f2,double f3,int nrec);

double adapt(
	double old_f(double),double a_old,double b_old,double acc,double eps);

double clenshaw_curtis(double old_f(double),double a_old,double b_old,double acc,double eps);
#endif
