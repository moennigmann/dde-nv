###############################################################################
#
# Npsol/ext_routines/debugging 
#
# revision history:
# 2008-02-28 Written  by dka 
#
###############################################################################

CCFLAGS= -g
CC= gcc

all: cwrap_npsol.o

clean: 
	rm -f *.o *.so *~* 

OBJS=funcJ0P.o funcF0P.o get_params.o extrhs.o extjac.o defs.o boundcond.o boundjac.o mcsubs.o opsubs.o f06subs.o cmsubs.o qrsubs.o lssolsubs.o npsolsubs.o chsubs.o srsubs.o set_params.o npsol_close.o  npsol_open.o npsolconfun.o npsolobjfun.o linalg_period.o zibconst.o zibsec.o period.o ftotides.o fortranjac.o fortranrhs.o twpbvpc.o callperiod.o infoprint.o	 

twpbvpc.o: twpbvpc.f
	${CC} -c ${CCFLAGS} twpbvpc.f -o twpbvpc.o

boundcond.o: boundcond.f
	${CC} -c ${CCFLAGS} boundcond.f -o boundcond.o

boundjac.o: boundjac.f
	${CC} -c ${CCFLAGS} boundjac.f -o boundjac.o
	
callperiod.o: callperiod.f
	${CC} -c ${CCFLAGS} callperiod.f -o callperiod.o
	
chsubs.o: chsubs.f
	${CC} -c ${CCFLAGS} chsubs.f -o chsubs.o	

defs.o: defs.f
	${CC} -c ${CCFLAGS} defs.f -o defs.o
	
extjac.o: extjac.f
	${CC} -c ${CCFLAGS} extjac.f -o extjac.o		

extrhs.o: extrhs.f
	${CC} -c ${CCFLAGS} extrhs.f -o extrhs.o	
	
f06subs.o: f06subs.f
	${CC} -c ${CCFLAGS} f06subs.f -o f06subs.o
	
cmsubs.o: cmsubs.f
	${CC} -c ${CCFLAGS} cmsubs.f -o cmsubs.o		

mcsubs.o: mcsubs.f
	${CC} -c ${CCFLAGS} mcsubs.f -o mcsubs.o

opsubs.o: opsubs.f
	${CC} -c ${CCFLAGS} opsubs.f -o opsubs.o

qrsubs.o: qrsubs.f
	${CC} -c ${CCFLAGS} qrsubs.f -o qrsubs.o

lssolsubs.o: lssolsubs.f
	${CC} -c ${CCFLAGS} lssolsubs.f -o lssolsubs.o

npsolsubs.o: lssolsubs.f
	${CC} -c ${CCFLAGS} npsolsubs.f -o npsolsubs.o

srsubs.o: srsubs.f
	${CC} -c ${CCFLAGS} srsubs.f -o srsubs.o

set_params.o: set_params.f
	${CC} -c ${CCFLAGS} set_params.f -o set_params.o
	
get_params.o: get_params.f
	${CC} -c ${CCFLAGS} get_params.f -o get_params.o	

npsol_close.o: npsol_close.f
	${CC} -c ${CCFLAGS} npsol_close.f -o npsol_close.o

npsol_open.o: npsol_open.f
	${CC} -c ${CCFLAGS} npsol_open.f -o npsol_open.o

npsolconfun.o: npsolconfun.f
	${CC} -c ${CCFLAGS} npsolconfun.f -o npsolconfun.o

npsolobjfun.o: npsolobjfun.f
	${CC} -c ${CCFLAGS} npsolobjfun.f -o npsolobjfun.o

period.o: period.f
	${CC} -c ${CCFLAGS} period.f -o period.o

linalg_period.o: linalg_period.f
	${CC} -c ${CCFLAGS} linalg_period.f -o linalg_period.o

zibconst.o: zibconst.f
	${CC} -c ${CCFLAGS} zibconst.f -o zibconst.o

zibsec.o: zibsec.f
	${CC} -c ${CCFLAGS} zibsec.f -o zibsec.o
	
ftotides.o: ftotides.f
	${CC} -c ${CCFLAGS} ftotides.f -o ftotides.o		

fortranjac.o: fortranjac.f
	${CC} -c ${CCFLAGS} fortranjac.f -o fortranjac.o
		
fortranrhs.o: fortranrhs.f
	${CC} -c ${CCFLAGS} fortranrhs.f -o fortranrhs.o
	
funcF0P.o: funcF0P.f
	${CC} -c ${CCFLAGS} funcF0P.f -o funcF0P.o
	
funcJ0P.o: funcJ0P.f
	${CC} -c ${CCFLAGS} funcJ0P.f -o funcJ0P.o

infoprint.o: infoprint.f
	${CC} -c ${CCFLAGS} infoprint.f -o infoprint.o				
	
cwrap_npsol.o: ${OBJS}
	${CC} -O2 ${OBJS} ${CCFLAGS} calltides.c scott.c cwrap_npsol.c -lTIDES -lm -lgfortran -o cwrap_npsol.o

