#include <stdio.h>

int main() {
	//Vars

    FILE *myFile;
    myFile = fopen("Rosenbrock.txt", "r");

    //read file into array
    int index[200];
    double x0[200];
    double x1[200];
    double f0[200];
    double f1[200];
    double f[200];
    int i=0;
    int status=1;
    printf(myFile);

    /*
    while (status != 0)
    {
        status=fscanf(myFile, "%d %d %d %d %d %n"\
        	, &index[i], &x0[i], &x1[i], &f0[i], &f1[i], &f[i]);
        printf("iter = %i \tx=%g \t y=%g\t f0 = %g\t f1 = %g\n"\
        	, index[i], x0[i], x1[i], f0[i], f1[i]);
        i++;
        

    }
	*/

    fclose(myFile);
return 0;
}