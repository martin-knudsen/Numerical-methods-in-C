CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main > out.txt

main: main.o Integrator.o

main.o Integrator.o: Integrator.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
