#ifndef HAVE_Monte_Carlo
#define HAVE_Monte_Carlo

void randomx(int dim,double *a,double *b,double *x);

void plainmc(int dim,double *a,double *b,
double f(double *x),int N,double *result,double *error);

#endif
