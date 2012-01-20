###############################################################################
#
# BVPSol/Makefile
#
# build from copy of NPSOL/Makefile
#
# target clean removes all *.o, *.so and *~* files from subdirectories to
#   BVPSol/ext_routines
#
# target all compiles all objects/shared objects in ext_routines/f2c
#
# revision history
#   2011-01-19 dka
# 
###############################################################################

clean:
	rm -f *~*
	make -C ./source clean
	make -C ./ext_routines/f2c/cwrap clean
	make -C ./ext_routines/f2c/shared_obj/period clean
	make -C ./ext_routines/f2c/shared_obj/twpbvp clean

all: 	
	make -C ./source all
	make -C ./ext_routines/f2c/cwrap all
	make -C ./ext_routines/f2c/shared_obj/period all
	make -C ./ext_routines/f2c/shared_obj/twpbvp all