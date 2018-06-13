linearshooting: linearshooting.o vec.o
	g++ linearshooting.o vec.o -o linearshooting 

linearshooting.o: linearshooting.cpp vec.h
	g++ -c linearshooting.cpp

vec.o: vec.cpp vec.h
	g++ -c vec.cpp

clean:
	rm *.o linearshooting 
