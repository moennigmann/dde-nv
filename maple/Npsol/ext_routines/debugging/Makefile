###############################################################################
#
# Npsol/ext_routines/debugging 
#
# revision history:
# 2008-02-28 Written  by dka 
#
###############################################################################

CCFLAGS= 
CC= gcc

all: check1.o

clean: 
	rm -f *.o *.so *~* 

OBJS=set_params.o npsolconfunad.o gnpsolconfun.o

set_params.o: set_params.f
	${CC} -c ${CCFLAGS} set_params.f -o set_params.o

npsolconfunad.o: npsolconfunad.c
	${CC} -c ${CCFLAGS} npsolconfunad.c -o npsolconfunad.o

gnpsolconfun.o: gnpsolconfun.f
	${CC} -c ${CCFLAGS} gnpsolconfun.f -o gnpsolconfun.o

check1.o: ${OBJS}
	${CC} ${OBJS} ${CCFLAGS} check1.c -lgfortran -o check1.o


