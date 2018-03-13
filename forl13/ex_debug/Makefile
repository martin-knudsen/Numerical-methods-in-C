CFLAGS = -w -std=gnu99
CFLAGS += 'gsl-config --cflags'
LDLIBS += $(gsl-config --libs)

all: output ; cat $@
output: main : $@ > $<
clean: ; $RM main output
