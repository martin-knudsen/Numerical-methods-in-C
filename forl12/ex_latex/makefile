CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: error_function
	
	./error_function -3 3 0.2 > error_function_data.txt

	gnuplot plot_error_function.gpi

error_function: error_function.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
