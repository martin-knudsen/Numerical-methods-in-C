CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 

all: Diagonalization
	./Diagonalization 

	gnuplot plot.gpi
Diagonalization: Diagonalization.o

.Phony: clean
clean: 
	$(RM) Diagonalization Diagonalization.o output.txt fit.txt 
