CFLAGS = -Wall -std=gnu99 `gsl-config --cflags` 
LDLIBS = `gsl-config --libs`  
 
all: QR_Gram_Schmidt

	./QR_Gram_Schmidt > output.txt

QR_Gram_Schmidt: QR_Gram_Schmidt.o

.Phony: clean 
clean: 
	$(RM) QR_Gram_Schmidt QR_Gram_Schmidt.o output.txt 
