CFLAGS = -Wall -std=gnu99
LDLIBS = -lm


.PHONEY: all
all: main-stdio main-cmdline 
	echo 1 2 3 4 5 | ./main-stdio > test.io.out.txt
	
	./main-cmdline 1 2 3 4 5 > test.cmd.out.txt

	echo `seq 0 0.02 6.5` > input.data.txt
	
	cat input.data.txt | ./main-stdio > out.io.txt
	./main-cmdline `cat input.data.txt` > out.cmd.txt

	gnuplot plot.gpi

main-stdio: main-stdio.o 

main-cmdline: main-cmdline.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete

