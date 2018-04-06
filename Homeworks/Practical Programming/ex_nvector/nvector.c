#include <stdio.h>
#include "nvector.h"
#include <assert.h>

nvector* nvector_alloc(int n){
  nvector* v = malloc(sizeof(nvector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
  return v;
}

void nvector_free(nvector* v){ free(v->data); free(v);}

void nvector_set(nvector* v, int i, double value){ 
	assert( 0 <= i && i < (*v).size );
	(*v).data[i]=value; 

}

double nvector_get(nvector* v, int i){return (*v).data[i]; }

double nvector_dot_product (nvector* u, nvector* v) {
	assert((*v).size==(*u).size);

	double result = 0;
	
	for(int i=0; i<(*u).size; i++) {
		
		result+= nvector_get (u, i)*nvector_get (v, i);
	}
	return result;
}

void nvector_print (char* s, nvector* v) {

	printf(s);
	printf(" [");
	for(int i=0; i<(*v).size; i++) {
		printf(" %g, ", (*v).data[i]);
	}
	printf("]\n");
}

void nvector_set_zero (nvector* v) {
	
	for(int i = 0; i<(*v).size; i++) {
		nvector_set(v, i, (double) 0.0);
	}
}

int nvector_equal (nvector* a, nvector* b) {
	assert((*a).size==(*b).size);
	for(int i = 0; i< (*a).size; i++) {
		if((*a).data[i]!=(*b).data[i]) {
			return 0;
		}
	}
	return 1;

}

void nvector_add (nvector* a, nvector* b) {
	if((*a).size!=(*b).size){printf("Vectors don't have same dimensions");}
	for(int i = 0; i< (*a).size; i++) {
		double value = (*a).data[i]+(*b).data[i];
		nvector_set(a, i, value);
	}
}

void nvector_sub (nvector* a, nvector* b) {
	if((*a).size!=(*b).size){printf("Vectors don't have same dimensions");}
	for(int i = 0; i< (*a).size; i++) {
		double value = (*a).data[i]-(*b).data[i];
		nvector_set(a, i, value);
	}
}

void nvector_scale (nvector* a, double x) {
	for(int i = 0; i < (*a).size; i++) {
		nvector_set(a, i, (*a).data[i]*x);
	}
}       


