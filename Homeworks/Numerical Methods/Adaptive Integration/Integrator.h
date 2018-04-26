#ifndef HAVE_Integrator
#define HAVE_Integrator

double adapt24(double f(double),double a,double b,
	double acc,double eps,double f2,double f3,int nrec);

double adapt(
	double f(double),double a,double b,double acc,double eps);

#endif
