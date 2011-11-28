#------------------------------------------------------------
#
# Aux/Makefile
#
# Top level makefile for Aux. 
#
# revision history:
# 2007-10-16 Written by M. Moennigmann on the way to the
#   ProcessNet 2007 conference. 
#
#------------------------------------------------------------
all: Aux.mpl

Aux.mpl: *.txt
	maple < Aux.txt 
	cat Aux.mpl AuxInitString.txt > tmp 
	mv -f tmp Aux.mpl 

clean:
	rm -f ./*~ ./BoxOperations/*~ ./CreateCurve/*~ ./FileOperations/*~ ./IntervalArithmetics/*~ ./LinearAlgebra/*~ ./LinearEqns/*~ ./ListOperations/*~ ./NLP/*~ ./Other/*~ ./Parsers/*~ ./Programming/*~ ./SystemClasses/*~ ./TransferNlpToGams/*~
	rm -f Aux.mpl
 