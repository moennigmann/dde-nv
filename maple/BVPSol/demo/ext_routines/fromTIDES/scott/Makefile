#############################################################################
#
# BVPSol/ext_routines/f2c/cwrap/scott/Makefile
#
#############################################################################

CCFLAGS= -c -shared -fPIC -O2
FFFLAGS= -c -shared -fPIC -O2
CC= gcc
FF= gfortran
LIBS=-lTIDES -lm

all: scott.o calltides.o
#
# compile c and fortran files
# 
scott.o: scott.c
	${CC} $(CCFLAGS) scott.c ${LIBS} -o scott.o

calltides.o: calltides.c
	${CC} $(CCFLAGS) calltides.c ${LIBS} -o calltides.o

clean: 
	rm -f *.o *~

