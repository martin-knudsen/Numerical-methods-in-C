CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 

all: linearSpline qspline cspline
	./linearSpline > output.txt
	./qspline 
	./cspline  
	gnuplot plot.gpi
	gnuplot plotcspline.gpi

cspline: cspline.o
qspline: qspline.o
linearSpline: linearSpline.o


.Phony: clean
clean: 
	$(RM) linearSpline linearSpline.o qspline qspline.o \
	cspline cspline.o cdata.txt csplinegsldata.txt \
	csplinetest.pdf lineardata.txt data_points.txt \
	output.txt qdata.txt splines.pdf 
