all:	
	doxygen Doxyfile
	cd latex && $(MAKE)
	mv latex/refman.pdf ./DDENLP-Documentation.pdf
	
	#rm *.m *.cpp
clean: 
	rm -r ./html/ ./latex/  *.m *.cpp
