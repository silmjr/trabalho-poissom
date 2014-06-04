CC=gcc
CFLAGS=-I./includes -lm -fopenmp -O3
OBJ=main.c libpoisson.o

# the $< is the first item in the dependencies list
# the $@ and $^ are the left and right sides
poisson: $(OBJ)
	$(CC) -o poisson $(OBJ) $(CFLAGS)

libpoisson.o: src/libpoisson.c
	$(CC) -c src/libpoisson.c $(CFLAGS)


clean:
	rm -f *.o badluk

	gcc -Wall -O3 -fopenmp -c main.c -o obj/Release/main.o
	gcc -Wall -O3 -fopenmp -c src/libpoisson.c -o obj/Release/src/libpoisson.o
	g++ -o bin/Release/poissonSerial obj/Release/main.o obj/Release/src/libpoisson.o -s -lgomp -lm

