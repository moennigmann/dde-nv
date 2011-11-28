###############################################################################
#
# Models/Makefile
#
###############################################################################

all: Models.mpl

Models.mpl: *.txt
	maple  < Models.txt 

clean:
	rm -f ./*~ Models.mpl
	make -C ./source clean
