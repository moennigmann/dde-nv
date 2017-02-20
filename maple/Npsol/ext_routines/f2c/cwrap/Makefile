###############################################################################
#
# NPSOL/ext_routines/f2c/cwrap/Makefile
# revision history:
# 2008-02-11 Makefile no longer copies files from ext_routines/cwrap,
# but all necessary files are stored in current directory. 
# 2008-02-11 Started with copy /media/disk/Backup/FirewireDisk-2007-04-08/
# backup/2006-07-backup-with-cp/home/mmo/maple/npsol - mmo
#
###############################################################################

CCFLAGS= -c -shared -O -fPIC 
FFFLAGS= -c -shared -O -fPIC
CC= gcc
FF= gfortran

all: npsol_open.o npsol_close.o cwrap_npsol.o npsolconfunad.o npsolobjfunad.o \
	set_params.o get_params.o

#
# compile c and fortran files
#
npsol_open.o: npsol_open.f
	${FF}  $(FFFLAGS) npsol_open.f -o npsol_open.o

npsol_close.o: npsol_close.f
	${FF}  $(FFFLAGS) npsol_close.f -o npsol_close.o

set_params.o: set_params.f
	${FF}  $(FFFLAGS) set_params.f -o set_params.o

get_params.o: get_params.f
	${FF}  $(FFFLAGS) get_params.f -o get_params.o

cwrap_npsol.o: cwrap_npsol.c
	${CC}  $(CCFLAGS) cwrap_npsol.c -o cwrap_npsol.o

npsolconfunad.o: npsolconfunad.c
	${CC}  $(CCFLAGS) npsolconfunad.c -o npsolconfunad.o 

npsolobjfunad.o: npsolobjfunad.c
	${CC}  $(CCFLAGS) npsolobjfunad.c -o npsolobjfunad.o 

clean: 
	rm -f *.o *~






