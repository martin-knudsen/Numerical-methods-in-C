#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


//create array of indices of largest elements in rows;

/* Built on top of the one taken from Dmitris chapter */
int jacobi_classic(gsl_matrix* A,gsl_vector* e,gsl_matrix* V, int* number_rot){
	/*Jacobi diagonalization;upper triangle of A is destroyed;
	e and V accumulate eigenvalues and eigenvectors*/

	// initiating the size, number of sweeps and if it has changed
	// in last sweep.
	int changed,sweeps=0,n=A->size1;

	// building a new array of largest indices ""
	gsl_vector* indices = gsl_vector_alloc(n);
	double new_el;
	double old_el;
	for(int i=0; i<n; i++){
		old_el = -INFINITY;		
		for(int j=i+1; j<n; j++){
			new_el = gsl_matrix_get(A, i, j);
			if(new_el>old_el) {
				old_el=new_el;
				gsl_vector_set(indices, i, j);
			}
		}
	}

	// put diagonal elements of A in vector e
	for(int i=0;i<n;i++)gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
		gsl_matrix_set_identity(V);
	
	// do while loop for performing sweeps.It will stop once there is no change 
	do{changed=0;sweeps++;int p,q;

	// iterating over all upper diagonal elements of A in a cyclic fashion
	for(p=0;p<n-1;p++){
		q = gsl_vector_get(indices, p);

		// Save all the old A values at p and q used in phi
		double app=gsl_vector_get(e,p);
		double aqq=gsl_vector_get(e,q);
		double apq=gsl_matrix_get(A,p,q);

		// save the angle that makes A_pq zero according to (3.10)
		double phi=0.5*atan2(2*apq,aqq-app);

		// save the cos and sin of phi as shortcuts
		double c=cos(phi),s=sin(phi);

		// calculate the new diagonal values at q and p of A1 according to (3.9)
		double app1=c*c*app-2*s*c*apq+s*s*aqq;
		double aqq1=s*s*app+2*s*c*apq+c*c*aqq;


		// If there is any change in the diagonal elements at q and p do
		if(app1!=app||aqq1!=aqq){changed=1; // mark that there was a change
			// update the new eigenvalues and set the appropriate zero'ed 
			// part of A, A_pq to zero
			gsl_vector_set(e,p,app1);
			gsl_vector_set(e,q,aqq1);
			gsl_matrix_set(A,p,q,0.0);
			*number_rot = *number_rot+1;
			// set the new value of p and q rows and columns at once
			// according to nr 2 and 3 from (3.9)
			for(int i=0;i<p;i++){
				double aip=gsl_matrix_get(A,i,p);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,i,p,c*aip-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*aip);}
			for(int i=p+1;i<q;i++){
				double api=gsl_matrix_get(A,p,i);
				double aiq=gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*api-s*aiq);
				gsl_matrix_set(A,i,q,c*aiq+s*api);}
			for(int i=q+1;i<n;i++){
				double api=gsl_matrix_get(A,p,i);
				double aqi=gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*api-s*aqi);
				gsl_matrix_set(A,q,i,c*aqi+s*api);}

			// It is worth noting that the first part of (3.9)
			// happens implicitly, because we use the same matrix 
			// and only the p, q rows and columns are affected by 
			// the rotation

			// Here we save the change on the eigenvectors
			// by the transformation according to (3.11)
			// again as before the V_ij are implicitly 
			// kept the same
			for(int i=0;i<n;i++){
				double vip=gsl_matrix_get(V,i,p);
				double viq=gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*vip-s*viq);
				gsl_matrix_set(V,i,q,c*viq+s*vip);}
	}}
	for(int i=0; i<n-1; i++){
		old_el = -1e8;		
		for(int j=i+1; j<n-1; j++){
			new_el = gsl_matrix_get(A, i, j);
			if(new_el>old_el) {
				old_el=new_el;
				gsl_vector_set(indices, i, j);
			}
		}
	}
	/*
	printm(A);
	printv(indices);
	*/
}
	while(changed!=0);
	return sweeps;}