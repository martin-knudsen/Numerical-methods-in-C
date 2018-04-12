#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
create array of indices of largest elements in rows;
do :
   for each row :
      eliminate the largest element of the row by Jacobi transformation;
      update the array of indices for the affected rows;
until converged.
*/