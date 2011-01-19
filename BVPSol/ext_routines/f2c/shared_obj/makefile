###############################################################################
#
# BVPSol/ext_routines/f2c/shared_obj/Makefile
#
# Prepares shared object from PERIOD code files used to link to maple. 
#
###############################################################################

CCFLAGS= -shared -fPIC
CC= gcc
FF= gfortran

all: period.so 

clean: 
	rm -f *.o *.so *~* 

OBJS= linalg_period.o zibconst.o zibsec.o period.o callperiod.o

period.o: period.f
	${FF} -c ${CCFLAGS} period.f -o period.o

linalg_period.o: linalg_period.f
	${FF} -c ${CCFLAGS} linalg_period.f -o linalg_period.o

zibconst.o: zibconst.f
	${FF} -c ${CCFLAGS} zibconst.f -o zibconst.o

zibsec.o: zibsec.f
	${FF} -c ${CCFLAGS} zibsec.f -o zibsec.o

callperiod.o: callperiod.f
	${FF} -c ${CCFLAGS} callperiod.f -o callperiod.o

period.so: $(OBJS) 
	$(CC) $(CCFLAGS) $(OBJS) -o period.so
