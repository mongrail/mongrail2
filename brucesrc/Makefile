CC = gcc
CFLAGS = -Wall -g

all: mongrail2

mongrail2: mongrail.o algorithms.o data.o inference.o error.o -lm
	$(CC) $(CFLAGS) -o mongrail2 mongrail.o algorithms.o data.o inference.o error.o -lm
mongrail.o: mongrail.c mongrail.h
	$(CC) $(CFLAGS) -c mongrail.c
algorithms.o: algorithms.c mongrail.h
	$(CC) $(CFLAGS) -c algorithms.c
data.o: data.c mongrail.h
	$(CC) $(CFLAGS) -c data.c
inference.o: inference.c mongrail.h
	$(CC) $(CFLAGS) -c inference.c
error.o: error.c mongrail.h
	$(CC) $(CFLAGS) -c error.c

clean:
	rm -f mongrail2 mongrail.o algorithms.o data.o inference.o error.o

