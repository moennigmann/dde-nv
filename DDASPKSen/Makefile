###############################################################################
#
# DDASPK/Makefile
#
# build from copy of NPSOL/Makefile
#
# target clean removes all *.o, *.so and *~* files from subdirectories to
#   NPSOL/ext_routines
#
# target all compiles all objects/shared objects in ext_routines/f2c
#
# target windows compiles all objects in ext_routines/f2c,
#    error messages in windows during executione of make can be ignored
#    (see remark in ext_routines/f2c/shared_obj/Makefile)
#
# revision history
#   080305 removed options solaris, win, common. Reduced to clean and all	    
#   050629 changed Macro MAPLE to maple.bat. This batch file must call the
#       desired maple command line version -jge
#   040727 added directory external_routines/win where now
#      compaq fortran compiler and MS c++ compiler are used 
#   032108 added target win and macro MAPLE which passes path to
#      maple executable -jge
#   010615 mmo
#
# To Do: 1) remove the make files in the win directory
#	 2) remove cleans for make -C ./ext_routines/cwrap clean.
# 
###############################################################################

clean:
	rm -f *~*
	make -C ./source clean
	#make -C ./ext_routines/cwrap clean
	#make -C ./ext_routines/shared_obj clean
	make -C ./ext_routines/f2c/cwrap clean
	make -C ./ext_routines/f2c/shared_obj clean

#	make -C ./ext_routines/win/cwrap clean
#	make -C ./ext_routines/win/fwrap clean
#	make -C ./ext_routines/win/shared_obj clean

all: 	
	make -C ./source all MAPLE=maple
	make -C ./ext_routines/f2c/cwrap all
	make -C ./ext_routines/f2c/shared_obj all