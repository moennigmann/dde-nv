all: 
	perl m2cpp.pl ../@VariableVector/VariableVector.m > ./VariableVector.cpp
	perl m2cpp.pl ../@DDE/DDE.m > ./DDE.cpp
	doxygen *.cpp
	rm *.cpp
clean:
	rm -r *.cpp ./html/ ./latex/
