CFLAGS=-fPIC 

all: cwrapper_test.so

cwrapper_test.o: cwrapper_test.c	
	gcc -c ${CFLAGS} cwrapper_test.c -o cwrapper_test.o

cwrapper_test.so: cwrapper_test.o
	gcc -shared ${CFLAGS} cwrapper_test.o -o cwrapper_test.so

clean:
	rm -f *.so *~ *.o






