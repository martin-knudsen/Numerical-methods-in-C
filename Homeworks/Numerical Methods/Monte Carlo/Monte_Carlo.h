#ifndef HAVE_Monte_Carlo
#define HAVE_Monte_Carlo

void randomx(int dim,double *a,double *b,double *x);

void plainmc(int dim,double *a,double *b,
double f(double *x),int N,double *result,double *error);

double adapt_2D(double old_f(double x, double y),double c(double), double d(double), double a_old,double b_old,double acc,double eps);

double adapt_3D(double old_f(double x, double y, double z),double ax,
 double bx, double ay,double by,double az, double bz,double acc,double eps);
#endif
