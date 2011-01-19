###############################################################################
#
# DDASPK/ext_routines/cwrap/Makefile
#
# revision history
#   011128 removed obsolete targets res.o, cwrap_ddaspk.so which were used 
#     for preliminary tests only; removed obsolete NPSOL comments
#   0110xx written by mkl using NPSOL/ext_routines/Makefile as a template
#
# Copyright (C) 2002 by M. Monnigmann for 
# Lehrstuhl f. Prozesstechnik,
# RWTH Aachen. All rights reserved.
#
###############################################################################

FFLAGS= -O -G
F77= f77 
CC= cc
LDFLAGS=

clean: 
	rm -f *.o *~*

all: cwrap_ddaspk.o jac_dummy.o psol_dummy.o g_res_dummy.o j_res_dummy.f

cwrap_ddaspk.o: cwrap_ddaspk.c
	${CC} -c $(FFLAGS) cwrap_ddaspk.c -o cwrap_ddaspk.o

jac_dummy.o: jac_dummy.f
	$(F77) -c jac_dummy.f -o jac_dummy.o

psol_dummy.o: psol_dummy.f
	$(F77) -c psol_dummy.f -o psol_dummy.o

g_res_dummy.o: g_res_dummy.f
	$(F77) -c g_res_dummy.f -o g_res_dummy.o


