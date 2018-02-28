CFLAGS = -Wall -std=gnu99 -I/usr/include

LDLIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm 

.PHONEY: all
all: logistic_diff equatorial_motion_i equatorial_motion_ii equatorial_motion_iii
	
	./logistic_diff
	./equatorial_motion_i
	./equatorial_motion_ii
	./equatorial_motion_iii

equatorial_motion_iii: equatorial_motion_iii.o	

equatorial_motion_ii: equatorial_motion_ii.o

equatorial_motion_i: equatorial_motion_i.o

logistic_diff: logistic_diff.o

.PHONEY: clean
clean:
	find ./ -type f -executable -delete
