CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	\time -f "%U" -ao time.txt ./main

	gnuplot plot.gpi

main: main.o jacobi_cyclic.o jacobi_eig_by_eig.o jacobi_classic.o QR.o QR_ls.o

main.o jacobi_cyclic.o jacobi_eig_by_eig.o jacobi_classic.o QR.o QR_ls.o: jacobi_cyclic.h jacobi_eig_by_eig.h jacobi_classic.h QR.h QR_ls.h

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
