#############################################################
#
# mapleDoc/Makefile
#
# Generates mapleDoc html files 
#
# revision history:
# 2011-11-23 Written by dka
#
#############################################################

all: ADIFOR AugSys2 Aux BVPSol Models Npsol BVPOpt
 
clean:
	rm -rf ./ADIFOR
	rm -rf ./AugSys2
	rm -rf ./Aux
	rm -rf ./BVPSol
	rm -rf ./Models
	rm -rf ./Npsol
	rm -rf ./BVPOpt

ADIFOR:
	java -jar ./mapledoc.jar ./../ADIFOR/source/ADIFOR.txt ./ADIFOR/ADIFOR.html

AugSys2:
	java -jar ./mapledoc.jar ./../AugSys2/AugSys2.txt ./AugSys2/AugSys2.html

Aux:
	java -jar ./mapledoc.jar ./../Aux/Aux.txt ./Aux/Aux.html

BVPSol:
	java -jar ./mapledoc.jar ./../BVPSol/source/BVPSol.txt ./BVPSol/BVPSol.html

Models:
	java -jar ./mapledoc.jar ./../Models/Models.txt ./Models/Models.html

Npsol:
	java -jar ./mapledoc.jar ./../Npsol/source/NPSOL.txt ./Npsol/Npsol.html

BVPOpt:
	java -jar ./mapledoc.jar ./../BVPOpt/BVPOpt.txt ./BVPOpt/BVPOpts.html 

