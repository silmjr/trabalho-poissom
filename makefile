
CC=gcc
CFLAGS=-I./includes -lm
OBJ=badluk.o flmoon.o julday.o nrutil.o

# the $< is the first item in the dependencies list
# the $@ and $^ are the left and right sides

%.o: %.c
        $(CC) -c -o $@ $< $(CFLAGS)

poisson: $(OBJ)
        $(CC) -o $@ $^ $(CFLAGS)

clean:
        rm -f *.o badluk

gcc -Wall -O3 -fopenmp -c main.c -o obj/Release/main.o
gcc -Wall -O3 -fopenmp -c src/libpoisson.c -o obj/Release/src/libpoisson.o
g++ -o bin/Release/poissonSerial obj/Release/main.o obj/Release/src/libpoisson.o -s -lgomp -lm

