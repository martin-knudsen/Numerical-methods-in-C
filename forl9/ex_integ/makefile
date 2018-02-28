CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: integral_example
	
	./integral_example

	
integral_example: integral_example.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
