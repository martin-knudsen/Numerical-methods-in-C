#ifndef HAVE_linearSpline
#define HAVE_linearSpline

int search1(double *x, double z);
double linterp(int n, double* x, double* y, double z);
double linterp_integ(int n, double *x, double *y, double z);
double linterp_der(int n, double *x, double *y, double z);
int main1();
#endif
