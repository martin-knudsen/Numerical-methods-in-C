CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main 

main: main.o QR_ls.o QR.o

main.o QR_ls.o QR.o: QR_ls.h QR.h

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
