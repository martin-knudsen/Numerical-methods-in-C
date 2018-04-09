CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 

all: linearSpline qspline cspline
	./linearSpline > output1.txt
	./qspline > output2.txt
	./cspline > output3.txt
	gnuplot plot.gpi

cspline: cspline.o
qspline: qspline.o
linearSpline: linearSpline.o


.Phony: clean
clean: 
	$(RM) linearSpline linearSpline.o \
	qspline qspline.o cspline cspline.o 
