CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main > out.txt
	 gnuplot plot.gpi

main: main.o Monte_Carlo.o Integrator.o QR.o QR_ls.o

main.o Monte_Carlo.o Integrator.o QR.o QR_ls.o: Monte_Carlo.h Integrator.h QR.h QR_ls.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
