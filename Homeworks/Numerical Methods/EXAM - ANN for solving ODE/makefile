CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main 
	 
	 

main: main.o ANN.o QR.o Newton.o 

main.o ANN.o QR.o Newton.o: ANN.h QR.h Newton.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
