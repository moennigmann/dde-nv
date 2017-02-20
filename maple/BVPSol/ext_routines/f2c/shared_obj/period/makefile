###############################################################################
#
# BVPSol/ext_routines/f2c/shared_obj/Makefile
#
# Prepares shared object from PERIOD code files used to link to maple. 
#
###############################################################################

CCFLAGS= -O2 -shared -fPIC
CC= gcc
FF= gfortran
LIBS=-lTIDES -lm

all: period.so 

clean: 
	rm -f *.o *.so *~* 

OBJS= linalg_period.o zibconst.o zibsec.o period.o ftotides.o callperiod.o

period.o: period.f
	${FF} -c ${CCFLAGS} period.f $(LIBS) -o period.o

linalg_period.o: linalg_period.f
	${FF} -c ${CCFLAGS} linalg_period.f $(LIBS) -o linalg_period.o

zibconst.o: zibconst.f
	${FF} -c ${CCFLAGS} zibconst.f $(LIBS) -o zibconst.o

zibsec.o: zibsec.f
	${FF} -c ${CCFLAGS} zibsec.f $(LIBS) -o zibsec.o

callperiod.o: callperiod.f
	${FF} -c ${CCFLAGS} callperiod.f $(LIBS) -o callperiod.o

ftotides.o: ftotides.f
	${FF} -c ${CCFLAGS} ftotides.f $(LIBS) -o ftotides.o	

period.so: $(OBJS) 
	$(CC) $(CCFLAGS) $(OBJS) -o period.so
