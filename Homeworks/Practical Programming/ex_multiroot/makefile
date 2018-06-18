CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: Rosenbrock shooting
	
	./Rosenbrock > out.txt
	./shooting
	
	gnuplot plot_hydrogen.gpi

Rosenbrock: Rosenbrock.o


shooting: shooting.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
