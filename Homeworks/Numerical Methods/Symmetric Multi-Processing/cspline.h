#ifndef HAVE_cspline
#define HAVE_cspline
typedef struct{int n; double *x,*y,*b,*c,*d;} cspline;
int search3(cspline* s, double z);
cspline* cspline_alloc(int n, double *x,double *y);
double cspline_eval(cspline *s,double z);
void cspline_free(cspline *s);
double cspline_derivative(cspline *s, double z);
double cspline_integral(cspline *s, double z);
int main3();
#endif
