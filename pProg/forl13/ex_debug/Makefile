CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`

all: main
	./main > output.txt
main: main.o
clean: ; $RM main output
