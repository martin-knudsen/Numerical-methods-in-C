#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "QR.h"
#define RND (double )rand()/RAND_MAX
#define FMT "%7.3f"

// These are all taken from Dmitri Fedorovs Numerical methods course
// from table 6.4-6.5
// performs reflection which "reflects" the highest point on the other side of the centroid
// p_re=p_ce+(p_ce-p_hi)
void reflection
(double *highest,double *centroid,int dim,double *reflected){
// for each dimension reflect the component of the highest point to the other side
	// of the corresponding centroid component
for(int i=0;i<dim;i++)reflected[i]=2*centroid[i]-highest[i];
}

// reflection of highest point and double the distance to the centroid
// p_ex=p_ce+2*(p_ce-p_hi)
void expansion
(double *highest,double *centroid,int dim,double *expanded){
// for all dimensions
for(int i=0;i<dim;i++)expanded[i]=3*centroid[i]-2*highest[i];
}

// highest point halves the distance to the centroid p_co=p_ce+1/2(p_hi-p_ce)
void contraction
(double *highest,double *centroid,int dim,double *contracted){
// for each dimension
for(int i=0;i<dim;i++)
contracted[i]=0.5*centroid[i]+0.5*highest[i];
}

// update all points to halve the distance to the lowest point
// p_k=1/2(p_k+p_lo)
void reduction(double**simplex,int dim,int lo){
// for all points
for(int k=0;k<dim+1;k++)
// if it is not the lowest point
if(k!=lo)
// halve the distance to the point dimension wise and put result in simplex
for(int i=0;i<dim;i++)
simplex[k][i]=0.5*(simplex[k][i]+simplex[lo][i]);
}

// just the euclidian distance between to vectors, th norm
double distance(double *a,double *b,int dim){
double s=0;for(int i=0;i<dim;i++)s+=pow(b[i]-a[i],2);
return sqrt(s);
}

// size of the simplex which is used to test if we have converged
// it measures the biggest distance between the points in the simplex
// if that distance is very small, it means that we are not changing much 
// and are converging. 
double size(double**simplex,int dim){
double s=0;for(int k=1;k<dim+1;k++){
double dist=distance(simplex[0],simplex[k],dim);
if(dist>s)s=dist;}
return s;
}

// this is just a method to find out what is the value place of 
// hi and lo in the simplex and sort the simplex according to the function values
void simplex_update(double**simplex,double *f_values,int d,
int *hi,int *lo,double *centroid){
*hi=0;*lo=0;double highest=f_values[0],lowest=f_values[0];
//find the integer of the highest and lowest point in simplex
for(int k=1;k<d+1;k++){
double next=f_values[k];
if(next>highest){highest=next;*hi=k;}
if(next<lowest){lowest=next;*lo=k;}}
// find the centroid of all points in simplex which are not the highest
// dimension wise. 
for(int i=0;i<d;i++){
double s=0;for(int k=0;k<d+1;k++)if(k!=*hi)s+=simplex[k][i];
centroid[i]=s/d;}
}

// gets the simplex and finds the function values of all points and updates
void simplex_initiate(
double fun(double *),double**simplex,double *f_values,int d,
int *hi,int *lo,double *centroid){
// finds the function values of all points
for(int k=0;k<d+1;k++)f_values[k]=fun(simplex[k]);
//find the new values of hi, lo, centroid
simplex_update(simplex,f_values,d,hi,lo,centroid);
}

/* This algorithm takes as input the function i question, a simplex,
which is basically just a collection of d+1 points in d input space
if our dimension is d. Here represented by a 2d array. Other input is
the dimension d, and the wished maximum simplex size. The algorithm works by reflection
expansion, contraction, reduction subroutines performed in a certain order.
The simplex is updated in the reduction step, which must be the "step".
The rest basically follows table 6.1.*/
int downhill_simplexx(
double F(double *),double**simplex,int d,double simplex_size_goal)
{
int hi,lo,k=0;double centroid[d],F_value[d+1],p1[d],p2[d];
// start by sorting and finding hi, lo, centroid, F_value
simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);
//convergence statement
while(size(simplex,d)>simplex_size_goal){
simplex_update(simplex,F_value,d,&hi,&lo,centroid);
//try reflection
reflection(simplex[hi],centroid,d,p1);double fre=F(p1);
// if reflection is good, do expansion
if(fre<F_value[lo]){//reflectionlooksgood:tryexpansion
expansion(simplex[hi],centroid,d,p2);double fex=F(p2);
if(fex<fre){//acceptexpansion
// apply expansion
for(int i=0;i<d;i++)simplex[hi][i]=p2[i];F_value[hi]=fex;}
else{//rejectexpansionandacceptreflection
// apply reflection
for(int i=0;i<d;i++)simplex[hi][i]=p1[i];F_value[hi]=fre;}}

else{//reflectionwasn’tgood
if(fre<F_value[hi]){//ok,acceptreflection
for(int i=0;i<d;i++)simplex[hi][i]=p1[i];F_value[hi]=fre;}
else{//trycontraction
contraction(simplex[hi],centroid,d,p1);double fco=F(p1);
if(fco<F_value[hi]){//acceptcontraction
for(int i=0;i<d;i++)simplex[hi][i]=p1[i];F_value[hi]=fco;}
// take step!
else{//doreduction
reduction(simplex,d,lo);
// again update and get F_values. 
simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);}}}
k++;}return k;
}