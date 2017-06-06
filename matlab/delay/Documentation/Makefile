all: 
	perl m2cpp.pl ../@VariableVector/VariableVector.m > ./VariableVector.cpp
	perl m2cpp.pl ../@DDE/DDE.m > ./DDE.cpp
	perl m2cpp.pl ../@EqualityConstraint/EqualityConstraint.m > ./EqualityConstraint.cpp
	perl m2cpp.pl ../@StStConstraint/StStConstraint.m > ./StStConstraint.cpp
	perl m2cpp.pl ../@NVConstraint/NVConstraint.m > ./NVConstraint.cpp
	doxygen *.cpp
	rm *.cpp
clean:
	rm -r ./html/ ./latex/
