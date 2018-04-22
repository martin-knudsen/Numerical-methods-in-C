CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main > out.txt

main: main.o Runge_Kutta.o 

main.o Runge_Kutta.o: Runge_Kutta.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
