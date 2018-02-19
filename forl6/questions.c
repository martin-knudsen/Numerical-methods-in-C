# include <stdio.h>

/* trace macro */
#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

int main() {
	return 0;
}