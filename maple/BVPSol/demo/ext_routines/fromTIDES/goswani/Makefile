#############################################################################
#
# BVPSol/ext_routines/f2c/cwrap/goswani/Makefile
#
#############################################################################

CCFLAGS= -c -shared -fPIC -O2
FFFLAGS= -c -shared -fPIC -O2
CC= gcc
FF= gfortran
LIBS=-lTIDES -lm

all: goswani.o calltides.o
#
# compile c and fortran files
# 
goswani.o: goswani.c
	${CC} $(CCFLAGS) goswani.c ${LIBS} -o goswani.o

calltides.o: calltides.c
	${CC} $(CCFLAGS) calltides.c ${LIBS} -o calltides.o

clean: 
	rm -f *.o *~

