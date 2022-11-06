CXX=g++ -std=c++17
OPT=-O3

PartiH2H: read.o index.o update.o main.o -lpthread -lboost_system -lboost_thread
	$(CXX) -g -o PartiCH read.o index.o update.o main.o -lpthread -lboost_system -lboost_thread

read.o:read.cpp
	$(CXX) -g -c $(OPT) read.cpp
index.o:index.cpp
	$(CXX) -g -c $(OPT) index.cpp
update.o:update.cpp
	$(CXX) -g -c $(OPT) update.cpp
main.o:main.cpp
	$(CXX) -g -c $(OPT) main.cpp

clean:
	rm *.o
	rm PartiCH

