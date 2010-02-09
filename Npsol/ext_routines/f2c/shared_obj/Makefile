###############################################################################
#
# Npsol/ext_routines/f2c/shared_obj/Makefile
#
# Prepares shared object from npsol code files used to link to maple. 
#
# revision history:
# 2008-02-11 Adjusted compiler flags. 
# 2008-02-11 Started with copy /media/disk/Backup/FirewireDisk-2007-04-08/
# backup/2006-07-backup-with-cp/home/mmo/maple/npsol - mmo
#
###############################################################################

CCFLAGS= -shared -fPIC
CC= gcc
FF= gfortran

all: npsol.so 

clean: 
	rm -f *.o *.so *~* 

OBJS=  blas1.o blas2.o mcsubs.o opsubs.o f06subs.o cmsubs.o qrsubs.o lssolsubs.o chsubs.o srsubs.o rtcomsubs.o nlssolsubs.o npcomsubs.o npsolsubs.o

blas1.o: blas1.f
	${FF} -c ${CCFLAGS} blas1.f -o blas1.o

blas2.o: blas2.f
	${FF} -c ${CCFLAGS} blas2.f -o blas2.o

mcsubs.o: mcsubs.f
	${FF} -c ${CCFLAGS} mcsubs.f -o mcsubs.o

opsubs.o: opsubs.f
	${FF} -c ${CCFLAGS} opsubs.f -o opsubs.o

f06subs.o: f06subs.f
	${FF} -c ${CCFLAGS} f06subs.f -o f06subs.o

cmsubs.o: cmsubs.f
	${FF} -c ${CCFLAGS} cmsubs.f -o cmsubs.o

qrsubs.o: qrsubs.f
	${FF} -c ${CCFLAGS} qrsubs.f -o qrsubs.o

lssolsubs.o: lssolsubs.f
	${FF} -c ${CCFLAGS} lssolsubs.f -o lssolsubs.o

chsubs.o: chsubs.f
	${FF} -c ${CCFLAGS} chsubs.f -o chsubs.o

srsubs.o: srsubs.f
	${FF} -c ${CCFLAGS} srsubs.f -o srsubs.o

rtcomsubs.o: rtcomsubs.f
	${FF} -c ${CCFLAGS} rtcomsubs.f -o rtcomsubs.o

nlssolsubs.o: nlssolsubs.f
	${FF} -c ${CCFLAGS} nlssolsubs.f -o nlssolsubs.o

npcomsubs.o: npcomsubs.f
	${FF} -c ${CCFLAGS} npcomsubs.f -o npcomsubs.o

npsolsubs.o: npsolsubs.f
	${FF} -c ${CCFLAGS} npsolsubs.f -o npsolsubs.o

npsol.so: $(OBJS) 
	$(CC) $(CCFLAGS) -fPIC $(OBJS) -o npsol.so

