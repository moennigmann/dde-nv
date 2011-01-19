###############################################################################
#
# DDASPK/ext_routines/f2c/shared_obj/Makefile
#
# copy of NPSOL/ext_routines/f2c/shared_obj/Makefile
#
# remark
#      error messages that occur in Windows while creating the shared
#      object file ddaspk.so because of the missing object
#      adf_dummy.o  in the linker command can be ignored, since maple
#      Procedure CreateSharedObject does not need ddaspk.so
#      if OperatingSystem=windows
# 
# revision history
#   080305 removed relicts and old comments. Changed compiler g77 against gfortran -tqu
#   030825 swapped order of ${CCFLAGS}$ ${OBJS}$ to ${OBJS}$ ${CCFLAGS}$ and
#   added -lf2c and -lm in CCFLAGS -jge 
#   020303 adf_dummy.f is no longer part of ddaspk.so but compiled into
#     separate object adf_dummy.o; adf_dummy.o is linked in to maple-ddaspk
#     interface by linker call in CreateSharedObject if option adifor is
#     not requested; 
#   011128 by mmo written to allow to use DDASPK on linux machine
#
# Copyright (C) 2002 by M. Monnigmann for 
# Lehrstuhl f. Prozesstechnik,
# RWTH Aachen. All rights reserved.
#
#ToDo: 1) remodel header
#
###############################################################################

GFORTRANFLAGS= -shared -O -fPIC 

DDASPK= adf_dummy.f  daux.f  ddaspk.f  dlinpk.f  dsensd.f  mpi_dummy.f

OBJS= daux.o  ddaspk.o  dlinpk.o  dsensd.o  mpi_dummy.o

clean: 
	rm -f *.c *.o *.so *~*

all: 
	gfortran $(GFORTRANFLAGS) $(DDASPK) -o ddaspk.so
	gfortran -c -o6 $(GFORTRANFLAGS) adf_dummy.f -o adf_dummy.o
