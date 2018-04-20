#ifndef HAVE_downhill_simplex
#define HAVE_downhill_simplex

void reflection
(double *highest,double *centroid,int dim,double *reflected);

void expansion
(double *highest,double *centroid,int dim,double *expanded);

void contraction
(double *highest,double *centroid,int dim,double *contracted);

void reduction(double **simplex,int dim,int lo);

double distance(double *a,double *b,int dim);

double size(double **simplex,int dim);

void simplex_update(double **simplex,double *f_values,int d,
int *hi,int *lo,double *centroid);

void simplex_initiate(
double fun(double *),double **simplex,double *f_values,int d,
int *hi,int *lo,double *centroid);

int downhill_simplex(
double F(double *),double**simplex,int d,double simplex_size_goal);
#endif
