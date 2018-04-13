CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: main

	 ./main >out.txt

main: main.o QR_ls.o

main.o QR_ls.o: QR_ls.h

.Phony: clean 
clean: 
	find . -type f -executable -delete
	rm -f *.o out*
