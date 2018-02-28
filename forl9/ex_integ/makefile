CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: integral_example energy_variational_principle
	
	./integral_example
	./energy_variational_principle > data.txt


energy_variational_principle             : energy_variational_principle.o normintegral.o hamiltonintegral.o        
energy_variational_principle.o normintegral.o hamiltonintegral.o: normintegral.h hamiltonintegral.h 


integral_example: integral_example.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
