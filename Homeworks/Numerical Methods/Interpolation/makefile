CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 

all: linearSpline qspline
	./linearSpline > output1.txt
	./qspline > output2.txt
	gnuplot plot.gpi

qspline: qspline.o
linearSpline: linearSpline.o

.Phony: clean
clean: 
	$(RM) linearSpline linearSpline.o  
