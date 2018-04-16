CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main 
	

main: main.o QR.o 

main.o QR.o: QR.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
