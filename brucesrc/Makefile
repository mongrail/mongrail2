CC = gcc
CFLAGS = -Wall -O3 -march=native -ffast-math
LDFLAGS = -lm -lhts

OBJS = mongrail.o algorithms.o data.o inference.o error.o vcf_reader.o

all: mongrail2

mongrail2: $(OBJS)
	$(CC) $(CFLAGS) -o mongrail2 $(OBJS) $(LDFLAGS)

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
vcf_reader.o: vcf_reader.c mongrail.h
	$(CC) $(CFLAGS) -c vcf_reader.c

clean:
	rm -f mongrail2 $(OBJS)

tidy:
	rm -f $(OBJS)
