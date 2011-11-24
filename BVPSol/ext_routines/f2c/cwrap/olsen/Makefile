###############################################################################
#
# BVPSol/ext_routines/f2c/cwrap/Makefile
#
###############################################################################

CCFLAGS= -c -shared -fPIC -O2
FFFLAGS= -c -shared -fPIC -O2
CC= gcc
FF= gfortran
LIBS=-lTIDES -lm

all: bvpopen.o bvpclose.o cwrap_bvpsol.o olsen.o calltides.o
#
# compile c and fortran files
#
bvpopen.o: bvpopen.f
	${FF}  $(FFFLAGS) bvpopen.f ${LIBS} -o  bvpopen.o

bvpclose.o: bvpclose.f
	${FF}  $(FFFLAGS) bvpclose.f ${LIBS} -o bvpclose.o

cwrap_bvpsol.o: cwrap_bvpsol.c
	${CC} $(CCFLAGS) cwrap_bvpsol.c ${LIBS} -o cwrap_bvpsol.o

olsen.o: olsen.c
	${CC} $(CCFLAGS) olsen.c ${LIBS} -o olsen.o

calltides.o: calltides.c
	${CC} $(CCFLAGS) calltides.c ${LIBS} -o calltides.o

clean: 
	rm -f *.o *~






