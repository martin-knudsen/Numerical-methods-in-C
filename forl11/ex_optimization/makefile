CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: rosenbrock_minimum
	
	./rosenbrock_minimum

	gnuplot rosenbrock_plot.gpi
	

rosenbrock_minimum: rosenbrock_minimum.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
