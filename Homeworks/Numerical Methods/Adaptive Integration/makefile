CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main 

main: main.o Integrator.o Runge_Kutta.o 

main.o Integrator.o Runge_Kutta.o: Integrator.h Runge_Kutta.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
