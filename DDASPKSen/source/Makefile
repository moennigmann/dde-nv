###############################################################################
#
# DDASPK/source/Makefile
#
# copy of NPSOL/source/Makefile
#
# Copyright (C) 2002 by M. Monnigmann for 
# Lehrstuhl f. Prozesstechnik,
# RWTH Aachen. All rights reserved.
#
# revision history
# 030821 added macro MAPLE that sets path to maple executable
###############################################################################

clean:
	rm -f DDASPK.mpl *~* Tests/*~*  SubsIntoSimpleWrapper/*~*
	rm -f SubsIntoSimpleWrapper/SimpleWrapper/*~*

all: DDASPK.mpl

DDASPK.mpl:
	maple < DDASPKSen.txt 
