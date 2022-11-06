CXX=g++ -std=c++17
OPT=-O3

partiPLL: Function.o Update.o PSL.o main.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o partiPLL Function.o Update.o PSL.o main.o -lpthread -lboost_system -lboost_thread

Function.o:Function.cpp
	$(CXX) -g -c $(OPT) Function.cpp
Update.o:Update.cpp
	$(CXX) -g -c $(OPT) Update.cpp
PSL.o:PSL.cpp
	$(CXX) -g -c $(OPT) PSL.cpp
main.o:main.cpp
	$(CXX) -g -c $(OPT) main.cpp

clean:
	rm *.o
	rm partiPLL

