###############################################################################
#
# ADIFOR/source/Makefile
#
# revision history
# 2006-08-15 switched order of targets to make all default target - mmo
#   030813 introduced macro MAPLE from which makefile receives path to
#   executable maple from calling makefile - jge  
#   020318 modelled after PITCON/source/Makefile
#
###############################################################################

all: ADIFOR.mpl

ADIFOR.mpl:
	maple < ADIFOR.txt

clean:
	rm -f ADIFOR.mpl *~

