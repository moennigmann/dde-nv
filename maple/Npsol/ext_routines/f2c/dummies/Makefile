###############################################################################
#
# revision history
#   010602 written after compilation of f77/c-routines failed on linux
#   because of linking problems
#
###############################################################################

FF= gfortran
CCFLAGS= -shared -O0 -fPIC 

clean:
	rm -f npsolconfun.o gnpsolconfun.o

all: npsolconfun.o 

npsolconfun.o: npsolconfun.f 
	${FF} -c ${CCFLAGS} npsolconfun.f -o npsolconfun.o
	${FF} -c ${CCFLAGS} gnpsolconfun.f -o gnpsolconfun.o


