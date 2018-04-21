#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* This method is partly taken from Dmitri Fedorovs lecture notes
but changed to perform 1 step instead of n steps.  
This is an implementation of the midpoint-Euler method step
It takes in the system of ODE's in the f - function, the current value 
of x, y(x) and the stepsize h as inputs. The output is the "integrated"
function y(x+h) and an error estimate dy. The ODE f gets input current x and y(x)
and outputs the righthand side of the ODE dydx=f(x,y)*/
void rkstep12(void f(double x,double* yx,double* dydx),
double x,double* yx,double h,double* yh,double* dy){
int i;double k0,yt,k12;
// calculate k0 using the ODE according to (7.14)
f(x,yx,k0);
// update using half-step temporary y-value
yt=yx+k0*h/2;
// calculate k1/2 using again (7.14)
f(x+h/2,yt,k12);
// update using (7.18) now k=k12. 
yh=yx+k12*h;
// estimate the error using (7.23) and looking at Butchers Tableau (7.25)
dy=(k0-k12)*h/2;
}

// This method is taken from Dmitri Fedorovs Numerical methods course lecture notes
int ode_driver(void f(int n,float x,float *y,float *dydx),
int n,float *xlist,float **ylist,
float b,float h,float acc,float eps,int max){
int i,k=0;float x,*y,s,err,normy,tol,a=xlist[0],yh[n],dy[n];
while(xlist[k]<b){
x=xlist[k],y=ylist[k];if(x+h>b)h=b-x;
rkstep12(f,n,x,y,h,yh,dy);
s=0;for(i=0;i<n;i++)s+=dy[i]*dy[i];err=sqrt(s);
s=0;for(i=0;i<n;i++)s+=yh[i]*yh[i];normy=sqrt(s);
tol=(normy*eps+acc)*sqrt(h/(b-a));
if(err<tol){/*accept step and continue*/
k++;if(k>max-1)return-k;/*uups*/
xlist[k]=x+h;for(i=0;i<n;i++)ylist[k][i]=yh[i];
}
if(err>0)h*=pow(tol/err,0.25)*0.95;elseh*=2;
}/*end while*/
return k+1;}/*return the number of entries in x list/ylist*/