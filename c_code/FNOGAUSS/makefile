CC = gcc
CFLAGS = -lm -L/usr/local/lib -llapack -lblas -lgsl -lgslcblas

calc : main.o diagonalize.o observable.o categorize.o
	$(CC) $(CFLAGS) main.o diagonalize.o observable.o categorize.o -o calc

%.o : %.c constants.h
	$(CC) $(CFLAGS) -c $< constants.h

clean :
	rm main.o diagonalize.o observable.o categorize.o
