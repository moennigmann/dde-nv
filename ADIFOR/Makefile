###############################################################################
#
# ADIFOR/Makefile
#
# revision history:
#   050629 changed Macro MAPLE to maple.bat. This batch file must call the
#          desired maple command line version 
#   030813 introduced target win and added Macro MAPLE which passes path to
#          maple executable to other Makefiles
#   020318 build from copy of PITCON/Makefile
#
###############################################################################

clean:
	rm -f *~*
	make -C ./source clean
common: 
	make -C ./source clean all MAPLE=maple

win:
	make -C ./source clean all MAPLE=maple.bat

Solaris: common

linux: common

all: common


