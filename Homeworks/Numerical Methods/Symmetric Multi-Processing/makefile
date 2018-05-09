CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 
CFLAGS += -fopenmp
LDLIBS += -lgomp
all: main

	 ./main >out.txt
	 gnuplot plot.gpi 
	 gnuplot plotcspline.gpi

main: main.o linearSpline.o qspline.o cspline.o  

main.o linearSpline.o qspline.o cspline.o: linearSpline.h qspline.h cspline.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
