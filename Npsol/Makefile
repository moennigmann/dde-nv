###############################################################################
#
# NPSOL/Makefile
#
# target clean removes all *.o, *.so and *~* files from subdirectories to
#   NPSOL/ext_routines
#
# target Solaris compiles all objects/shared object in ext_routines
#   except for subdirs of ext_routines/f2c which are used under linux only
#
# target linux compiles all objects/shared objects in ext_routines/f2c
#
# revision history
# 030418 removed dir debugging from clean target
# 021217 added directory debugging to clean target
# 010615 written by mmo
#
###############################################################################

clean:
	rm -f *~
	make -C ./source clean
	make -C ./ext_routines/cwrap clean
	make -C ./ext_routines/shared_obj clean
	make -C ./ext_routines/dummies clean
	make -C ./ext_routines/f2c/cwrap clean
	make -C ./ext_routines/f2c/shared_obj clean
	make -C ./ext_routines/f2c/dummies clean

common: 
	make -C ./source all

Solaris: common
	make -C ./ext_routines/cwrap all
	make -C ./ext_routines/shared_obj all
	make -C ./ext_routines/dummies all

linux: common
	make -C ./ext_routines/f2c/cwrap all
	make -C ./ext_routines/f2c/shared_obj all
	make -C ./ext_routines/f2c/dummies all

