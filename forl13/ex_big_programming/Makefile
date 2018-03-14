CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs` 

all: Diagonalization
	./Diagonalization > output.txt

	gnuplot plot.gpi
Diagonalization: Diagonalization.o
clean: ; $RM main output
