CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: integral_example energy_variational_principle the_reason
	
	./integral_example
	./energy_variational_principle > data.txt

	gnuplot plot.gpi

	./the_reason



energy_variational_principle             : energy_variational_principle.o normintegral.o hamiltonintegral.o        
energy_variational_principle.o normintegral.o hamiltonintegral.o: normintegral.h hamiltonintegral.h 


integral_example: integral_example.o

the_reason: the_reason.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
