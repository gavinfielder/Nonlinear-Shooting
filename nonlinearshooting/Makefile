nonlinearshooting: nonlinearshooting.o shootingIterator.o
	g++ nonlinearshooting.o shootingIterator.o -o nonlinearshooting 

nonlinearshooting.o: nonlinearshooting.cpp shootingIterator.h
	g++ -c nonlinearshooting.cpp

shootingIterator.o: shootingIterator.cpp shootingIterator.h
	g++ -c shootingIterator.cpp

clean:
	rm *.o nonlinearshooting 
