#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

/* These methods are inspired Dmitri Fedorovs lecture notes but implemented using gsl_vectors
and gsl_matrices and gsl_blas. 
This is an implementation of the midpoint-Euler method step
It takes in the system of ODE's in the f - function, the current value 
of x, y(x),the stepsize h and the dimension of the problem n as inputs. The output is the "integrated"
function y(x+h) and an error estimate dy. The ODE f gets input current x and y(x)
and outputs the righthand side of the ODE dydx=f(x,y)*/
void rkstep12(void f(int n,double x,gsl_vector *y,gsl_vector *dydx),
int n, double x,gsl_vector* yx,double h,gsl_vector* yh,gsl_vector* dy){
int i;
gsl_vector* k0=gsl_vector_alloc(n);
gsl_vector* yt=gsl_vector_alloc(n);
gsl_vector* k12=gsl_vector_alloc(n);
// calculate k0 using the ODE according to (7.14)
f(n,x,yx,k0);
// update using half-step temporary y-value for each equation
for(i=0;i<n;i++) {
	gsl_vector_set(yt,i,gsl_vector_get(yx,i)+gsl_vector_get(k0,i)*h/2);
}
// calculate k1/2 using again (7.14)
f(n,x+h/2,yt,k12);
// update using (7.18) now k=k12 for each equation
for(i=0;i<n;i++) {
	gsl_vector_set(yh,i,gsl_vector_get(yx,i)+gsl_vector_get(k12,i)*h);
}
// estimate the error using (7.23) and looking at Butchers Tableau (7.25)
for(i=0;i<n;i++) 
	gsl_vector_set(dy,i,(gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*h/2);
gsl_vector_free(k0);
gsl_vector_free(yt);
gsl_vector_free(k12);
}

/* Input is the ODE system as described above, the dimension of the problem n, 
the endpoint b, stepsize h, absolute accuracy acc, rel. accuracy eps, 
max is the maximum number of allowed steps. Output is the number of steps taken,
k, the recorded path in the form of x values and y(x) values of all steps taken
in the form of a 1d array xlist and a 2d-array ylist*/
int ode_driver(void f(int n,double x,gsl_vector *y,gsl_vector *dydx),
int n,gsl_vector *xlist,gsl_matrix *ylist,
double b,double h,double acc,double eps,int max){
/*k is the current stepnumber, x the current x value, y an array of the current
y-values, err,normy,tol the current error, |y|, tolerance, a the startpoint
of the iteration, yh, dy 1d-arrays of current yh and dy values at step. */
int i,k=0;double x,err,normy,tol,a=gsl_vector_get(xlist,0);
gsl_vector* y=gsl_vector_alloc(n);
gsl_vector* yh=gsl_vector_alloc(n);
gsl_vector* dy=gsl_vector_alloc(n);
//keep taking steps until value of x is larger than wanted endvalue b
while(gsl_vector_get(xlist,k)<b){
//current x value and y values.   
x=gsl_vector_get(xlist,k); gsl_matrix_get_row(y,ylist,k);
//if the next step steps over the endpoint b, choose step size that exactly make it to b 
if(x+h>b)h=b-x;
//take the step here with rk12 algorithm
rkstep12(f,n,x,y,h,yh,dy);
//find the new estimate of error and euclidian norm(y) of the step according to (7.41). 
err=gsl_blas_dnrm2(dy);
normy=gsl_blas_dnrm2(yh);
// also use (7.41) to calculate the accepted tolerance

tol=(normy*eps+acc)*sqrt(h/(b-a));
// if the error is within the acceptable tolerance accept step
if(err<tol){/*accept step and continue*/
k++;
// if we have taken the maximum number of steps return negative number, which signales we didn't find
// the solution within the relative and absolute accuracy and maximum steps.
if(k>max-1)return-k;/*uups*/
// save the step to the paths
gsl_vector_set(xlist,k,x+h);
for(i=0;i<n;i++) gsl_matrix_set(ylist,k,i,gsl_vector_get(yh,i));
}
// there is an error, update according to (7.40). if there is no error 
// just take the same step length!
if(err>0) h*=pow(tol/err,0.25)*0.95; else h*=2;
}/*end while*/
return k;

gsl_vector_free(y);
gsl_vector_free(yh);
gsl_vector_free(dy);


}/*return the number of entries in x list/ylist*/