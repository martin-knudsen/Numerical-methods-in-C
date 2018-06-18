#include <stdio.h>

int main() {
	FILE *myfile;
myfile = fopen("out.txt", "a+");
	fprintf(myfile,"The reason alfa=0.5 is a minimum is that this is the theoretical value for the QM harmonic oscillator.\n");
	fclose(myfile);	
	return 0;
}
