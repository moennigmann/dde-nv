#############################################################################
#
# BVPSol/ext_routines/f2c/cwrap/olsen/Makefile
#
#############################################################################

CCFLAGS= -c -shared -fPIC -O2
FFFLAGS= -c -shared -fPIC -O2
CC= gcc
FF= gfortran
LIBS=-lTIDES -lm

all: olsen.o calltides.o
#
# compile c and fortran files
# 
olsen.o: olsen.c
	${CC} $(CCFLAGS) olsen.c ${LIBS} -o olsen.o

calltides.o: calltides.c
	${CC} $(CCFLAGS) calltides.c ${LIBS} -o calltides.o

clean: 
	rm -f *.o *~

