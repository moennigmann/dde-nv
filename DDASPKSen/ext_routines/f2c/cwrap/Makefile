###############################################################################
#
# DDASPK/ext_routines/f2c/cwrap/Makefile
#
# revision history
# 080305 removed flags -M64; they are set on default; changed fortran 
#	 compiler g77 against gfortran 
# 011128 copy of NPSOL/ext_routines/f2c/cwrap/Makefile;
#   subroutine names may not contain "_", sed commands below remove
#   underline characters in all source files, replacements are collected
#   in macro SEDARGS below
#
# Copyright (C) 2002 by M. Monnigmann for 
# Lehrstuhl f. Prozesstechnik,
# RWTH Aachen. All rights reserved.
#
###############################################################################

CCFLAGS= -c -shared -O -fPIC
GFORTRANFLAGS= -c -shared -O -fPIC 
CC= gcc

clean: 
	rm -f *.o *~* *.f *.c 

all: cwrap_ddaspk.o psol_dummy.o g_res_dummy.o jac_dummy.o #j_res_dummy.o 

#
# compile c files
#

psol_dummy.o:
	cp -f ../../cwrap/psol_dummy.f .
	gfortran $(GFORTRANFLAGS) psol_dummy.f -o psol_dummy.o

jac_dummy.o: 
	cp -f ../../cwrap/jac_dummy.f .
	gfortran $(GFORTRANFLAGS) jac_dummy.f -o jac_dummy.o

g_res_dummy.o: 
	cp -f ../../cwrap/g_res_dummy.f .
	gfortran $(GFORTRANFLAGS) g_res_dummy.f -o g_res_dummy.o

cwrap_ddaspk.o: 
	cp -f ../../cwrap/cwrap_ddaspk.c .
	${CC} $(CCFLAGS) cwrap_ddaspk.c -o cwrap_ddaspk.o

