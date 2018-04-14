CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main > out.txt
	 gnuplot plot.gpi

main: main.o QR_ls.o QR.o s_ls.o jacobi_cyclic.o

main.o QR_ls.o QR.o s_ls.o jacobi_cyclic.o: QR_ls.h QR.h s_ls.h jacobi_cyclic.h

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
