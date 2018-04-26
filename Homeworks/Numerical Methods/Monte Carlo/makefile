CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main 

main: main.o Monte_Carlo.o Integrator.o 

main.o Monte_Carlo.o Integrator.o: Monte_Carlo.h Integrator.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
