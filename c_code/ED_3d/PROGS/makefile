CC=icc
LIBS=-lgsl -lgslcblas -lm
CFLAGS= -std=c99 -Wall
OBJS=u1link.o construct.o common.o hamiltonian.o hamiltonianF.o

all: link3du1

link3du1: ${OBJS}
	${CC} -o link3du1 ${OBJS} ${LIBS}

clean:
	rm -rf link3du1 ${OBJS} *.optrt
