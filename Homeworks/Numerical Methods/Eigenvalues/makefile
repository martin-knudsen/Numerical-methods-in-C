CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	\time -f "%U" -ao time.txt ./main > output.txt

main: main.o jacobi_cyclic.o jacobi_eig_by_eig.o jacobi_classic.o
	
main.o jacobi_cyclic.o:jacobi_cyclic.h
main.o jacobi_eig_by_eig.o:jacobi_eig_by_eig.h
main.o jacobi_classic.o:jacobi_classic.h

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
