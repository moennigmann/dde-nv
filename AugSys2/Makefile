#------------------------------------------------------------
#
# AugSys2/Makefile
#
# Top level makefile for AugSys2. 
#
# revision history:
# 2007-12-04 Written by dka
#
#------------------------------------------------------------
all: AugSys2.mpl

AugSys2.mpl: *.txt
	maple  < AugSys2.txt 

clean:
	rm -f ./*~ ./NeimarkSackerNV/*~ ./FlipAndFoldNV/*~ ./ConstraintInitializer/*~
	rm -f ./ConstraintInitializer/CriticalPoint*~ ./ConstraintInitializer/FixPoint*~
	rm -f AugSys2.mpl

