###############################################################################
#
# BVPSol/ext_routines/f2c/cwrap/Makefile
#
###############################################################################

CCFLAGS= -c -shared -O -fPIC 
FFFLAGS= -c -shared -O -fPIC
CC= gcc
FF= gfortran

all: bvpopen.o bvpclose.o cwrap_bvpsol.o
#
# compile c and fortran files
#
bvpopen.o: bvpopen.f
	${FF}  $(FFFLAGS) bvpopen.f -o bvpopen.o

bvpclose.o: bvpclose.f
	${FF}  $(FFFLAGS) bvpclose.f -o bvpclose.o

cwrap_bvpsol.o: cwrap_bvpsol.c
	${CC}  $(CCFLAGS) cwrap_bvpsol.c -o cwrap_bvpsol.o

clean: 
	rm -f *.o *~






