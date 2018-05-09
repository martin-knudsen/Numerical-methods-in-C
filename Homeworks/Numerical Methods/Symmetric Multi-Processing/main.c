#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "linearSpline.h"
#include "qspline.h"
#include "cspline.h"

int main(void) {
	
	printf("Multiprocessing of Interpolation exercise.\n"
		"This is made by parallel processing the linear, quadric and cubic spline\n"
		"part of that exercise. The results are the same but each output is in \n"
		"out1.txt, out2.txt and out3.txt respectively");
		
	// This template is inspired by Dmitris implementation of openMP
	// it works by running the 3 main functions in parallel
	#pragma omp parallel sections
		//the following sections wil be run parallelly in separate threads
		{
		#pragma omp section //first thread will run this block of code
			{
			main1();
			}
		#pragma omp section //second thread will run this block of code
			{
			main2();
			}
		#pragma omp section //third thread will run this block of code
			{
			main3();
			}
		}

	return EXIT_SUCCESS;
}