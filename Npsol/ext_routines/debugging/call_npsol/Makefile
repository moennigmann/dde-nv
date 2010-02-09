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

all: cwrap_npsol.o

clean: 
	rm -f *.o *.so *~* 

OBJS=mcsubs.o opsubs.o f06subs.o cmsubs.o qrsubs.o lssolsubs.o npsolsubs.o chsubs.o srsubs.o set_params.o npsol_close.o  npsol_open.o npsolconfun.o npsolobjfun.o

mcsubs.o: mcsubs.f
	${CC} -c ${CCFLAGS} mcsubs.f -o mcsubs.o

opsubs.o: opsubs.f
	${CC} -c ${CCFLAGS} opsubs.f -o opsubs.o

f06subs.o: f06subs.f
	${CC} -c ${CCFLAGS} f06subs.f -o f06subs.o

cmsubs.o: cmsubs.f
	${CC} -c ${CCFLAGS} cmsubs.f -o cmsubs.o

qrsubs.o: qrsubs.f
	${CC} -c ${CCFLAGS} qrsubs.f -o qrsubs.o

lssolsubs.o: lssolsubs.f
	${CC} -c ${CCFLAGS} lssolsubs.f -o lssolsubs.o

npsolsubs.o: lssolsubs.f
	${CC} -c ${CCFLAGS} npsolsubs.f -o npsolsubs.o

chsubs.o: chsubs.f
	${CC} -c ${CCFLAGS} chsubs.f -o chsubs.o

srsubs.o: srsubs.f
	${CC} -c ${CCFLAGS} srsubs.f -o srsubs.o

set_params.o: set_params.f
	${CC} -c ${CCFLAGS} set_params.f -o set_params.o

npsol_close.o: npsol_close.f
	${CC} -c ${CCFLAGS} npsol_close.f -o npsol_close.o

npsol_open.o: npsol_open.f
	${CC} -c ${CCFLAGS} npsol_open.f -o npsol_open.o

npsolconfun.o: npsolconfun.f
	${CC} -c ${CCFLAGS} npsolconfun.f -o npsolconfun.o

npsolobjfun.o: npsolobjfun.f
	${CC} -c ${CCFLAGS} npsolobjfun.f -o npsolobjfun.o

cwrap_npsol.o: ${OBJS}
	${CC} ${OBJS} ${CCFLAGS} cwrap_npsol.c -lgfortran -o cwrap_npsol.o

