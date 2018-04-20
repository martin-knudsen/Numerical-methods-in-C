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
void reflection
(double *highest,double *centroid,int dim,double *reflected){
for(int i=0;i<dim;i++)reflected[i]=2*centroid[i]-highest[i];
}
void expansion
(double *highest,double *centroid,int dim,double *expanded){
for(int i=0;i<dim;i++)expanded[i]=3*centroid[i]-2*highest[i];
}
void contraction
(double *highest,double *centroid,int dim,double *contracted){
for(int i=0;i<dim;i++)
contracted[i]=0.5*centroid[i]+0.5*highest[i];
}
void reduction(double**simplex,int dim,int lo){
for(int k=0;k<dim+1;k++)if(k!=lo)for(int i=0;i<dim;i++)
simplex[k][i]=0.5*(simplex[k][i]+simplex[lo][i]);
}
double distance(double *a,double *b,int dim){
double s=0;for(int i=0;i<dim;i++)s+=pow(b[i]-a[i],2);
return sqrt(s);
}
double size(double**simplex,int dim){
double s=0;for(int k=1;k<dim+1;k++){
double dist=distance(simplex[0],simplex[k],dim);
if(dist>s)s=dist;}
return s;
}

void simplex_update(double**simplex,double *f_values,int d,
int *hi,int *lo,double *centroid){
*hi=0;*lo=0;double highest=f_values[0],lowest=f_values[0];
for(int k=1;k<d+1;k++){
double next=f_values[k];
if(next>highest){highest=next;*hi=k;}
if(next<lowest){lowest=next;*lo=k;}}
for(int i=0;i<d;i++){
double s=0;for(int k=0;k<d+1;k++)if(k!=*hi)s+=simplex[k][i];
centroid[i]=s/d;}
}
void simplex_initiate(
double fun(double *),double**simplex,double *f_values,int d,
int *hi,int *lo,double *centroid){
for(int k=0;k<d+1;k++)f_values[k]=fun(simplex[k]);
simplex_update(simplex,f_values,d,hi,lo,centroid);
}

int downhill_simplex(
double F(double *),double**simplex,int d,double simplex_size_goal)
{
int hi,lo,k=0;double centroid[d],F_value[d+1],p1[d],p2[d];
simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);
while(size(simplex,d)>simplex_size_goal){
simplex_update(simplex,F_value,d,&hi,&lo,centroid);
reflection(simplex[hi],centroid,d,p1);double fre=F(p1);
if(fre<F_value[lo]){//reflectionlooksgood:tryexpansion
expansion(simplex[hi],centroid,d,p2);double fex=F(p2);
if(fex<fre){//acceptexpansion
for(int i=0;i<d;i++)simplex[hi][i]=p2[i];F_value[hi]=fex;}
else{//rejectexpansionandacceptreflection
for(int i=0;i<d;i++)simplex[hi][i]=p1[i];F_value[hi]=fre;}}
else{//reflectionwasn’tgood
if(fre<F_value[hi]){//ok,acceptreflection
for(int i=0;i<d;i++)simplex[hi][i]=p1[i];F_value[hi]=fre;}
else{//trycontraction
contraction(simplex[hi],centroid,d,p1);double fco=F(p1);
if(fco<F_value[hi]){//acceptcontraction
for(int i=0;i<d;i++)simplex[hi][i]=p1[i];F_value[hi]=fco;}
else{//doreduction
reduction(simplex,d,lo);
simplex_initiate(F,simplex,F_value,d,&hi,&lo,centroid);}}}
k++;}return k;
}