CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main > out.txt

	 gnuplot plot.gpi
	

main: main.o Newton.o QR.o downhill_simplex.o

main.o Newton.o QR.o downhill_simplex.o: Newton.h QR.h downhill_simplex.h

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
