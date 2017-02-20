###############################################################################
#
# BVPsol/ext_subroutines/debugging 
#
# revision history:
# 2011-01-14 written  by dka 
#
###############################################################################

CCFLAGS= -g
CC= gcc

all: cwrap_bvpsol.o 

clean: 
	rm -f *.o *.so *~* 

OBJS=linalg_period.o zibconst.o zibsec.o fortranrhs.o fortranjac.o period.o callperiod.o

zibconst.o: zibconst.f
	${CC} -c ${CCFLAGS} zibconst.f -o zibconst.o

zibsec.o: zibsec.f
	${CC} -c ${CCFLAGS} zibsec.f -o zibsec.o

fortranrhs.o: fortranrhs.f
	${CC} -c ${CCFLAGS} fortranrhs.f -o fortranrhs.o

fortranjac.o: fortranjac.f
	${CC} -c ${CCFLAGS} fortranjac.f -o fortranjac.o

period.o: period.f
	${CC} -c ${CCFLAGS} period.f -o period.o

callperiod.o: callperiod.f
	${CC} -c ${CCFLAGS} callperiod.f -o callperiod.o

linalg_period.o: linalg_period.f
	${CC} -c ${CCFLAGS} linalg_period.f -o linalg_period.o

cwrap_bvpsol.o: ${OBJS}
	${CC} ${OBJS} ${CCFLAGS} cwrap_bvpsol.c -lgfortran -o cwrap_bvpsol.o

