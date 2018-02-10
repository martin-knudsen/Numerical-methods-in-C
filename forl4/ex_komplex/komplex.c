#include "komplex.h"
#include <stdio.h>
#include <float.h>
#include <math.h>

void komplex_print (char *s, komplex a) {
	printf ("%s (%g,%g)\n", s, a.re, a.im);
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;
}

int komplex_equal (komplex a, komplex b) {
	if(a.re==b.re && a.im==b.im){
		return 1;
	}
	else {
		return 0;
	}
}
	
komplex komplex_mul (komplex a, komplex b) {
	komplex result = {a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re};
	return result;
}

komplex komplex_conjugate(komplex z) {
	komplex result = {z.re, -z.im};
	return result;

}

komplex komplex_div (komplex a, komplex b) {
	komplex b_conj = komplex_conjugate(b);
	komplex num = komplex_mul(a, b_conj);
	komplex denom = komplex_mul(b, b_conj);
	komplex result = {num.re/denom.re, num.im/denom.re};
	return result;
}


komplex komplex_abs (komplex z) {
	komplex result = {pow(z.re,2)+pow(z.im,2), 0};
	return result;
}

komplex komplex_exp (komplex z) {
	komplex exp_im_part = {cos(z.im), sin(z.im)};
	komplex exp_re_part = {exp(z.re) , 0};
	komplex result = komplex_mul(exp_re_part,exp_im_part);
	return result;  
}

komplex komplex_sin (komplex z) {
	komplex iz = {-z.im, z.re};
	komplex exp1 = komplex_exp(iz);

	komplex iz_neg = {z.im, -z.re};
	komplex exp2 = komplex_exp(iz_neg);
	komplex neg_term ={-exp2.re,exp2.im};  
	
	komplex num = komplex_add(exp1,neg_term);
	komplex denom = {0, 2};

	komplex result = komplex_div(num, denom);
	return result;
}

/* ... */