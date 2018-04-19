CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main >out.txt
	

main: main.o Newton.o QR.o 

main.o Newton.o QR.o: Newton.h QR.h 

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
