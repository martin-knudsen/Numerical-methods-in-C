CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm

.PHONEY: all
all: Airy_stdio 
	

	echo `seq -10 0.02 2` > input.data.txt
	
	cat input.data.txt | ./Airy_stdio > out.io.txt

	gnuplot plot.gpi

Airy_stdio: Airy_stdio.o


.PHONEY: clean
clean:
	find ./ -type f -executable -delete
