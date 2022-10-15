CC=mpicc
CFLAGS=-I. -W -pg -O3 -ffast-math -funroll-loops
DEPS=ising.h
OBJ=ising.o isingfunc.o
LIBS=-lm

%.o:%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ising.x: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	
.PHONY:clean

clean:
	rm -f *.o *.x *.out *.dat *~ core
