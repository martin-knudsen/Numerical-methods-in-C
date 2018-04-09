CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 

all: linearSpline
	./linearSpline > output.txt

	gnuplot plot.gpi

linearSpline: linearSpline.o

.Phony: clean
clean: 
	$(RM) linearSpline linearSpline.o  
