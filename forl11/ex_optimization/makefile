CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: rosenbrock_minimum least_squares
	
	./rosenbrock_minimum

	gnuplot rosenbrock_plot.gpi
	
	./least_squares

least_squares: least_squares.o
rosenbrock_minimum: rosenbrock_minimum.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
