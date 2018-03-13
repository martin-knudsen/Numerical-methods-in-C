CFLAGS = -Wall -std=c99 
LDLIBS = -lm

.PHONEY: all clean indent

all              : output.txt          ; cat $<           # "$<" is the first prerequisite
output.txt       : main             ; ./$< > $@        # "$@" is the target
main             : main.o nvector.o # built-in linking rule is used here         
main.o nvector.o : nvector.h        # built-in compilation rule is used here

clean            :                  ; $(RM) main out* *.o *~
