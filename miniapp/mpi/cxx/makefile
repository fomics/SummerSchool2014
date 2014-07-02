CXX=CC
ifeq ($(PE_ENV),GNU)
    # the gnu compiler requires that we specifically enable OpenMP
	#CXXFLAGS=-O0 -g -fopenmp -DDEBUG
	CXXFLAGS=-O3 -fopenmp
endif
ifeq ($(PE_ENV),CRAY)
	CXXFLAGS=-O3
endif

SOURCES = stats.cpp data.cpp operators.cpp linalg.cpp
HEADERS = stats.h   data.h   operators.h   linalg.h
OBJ     = stats.o   data.o   operators.o   linalg.o

.SUFFIXES: .cpp

all: main

stats.o: stats.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c stats.cpp -o stats.o

data.o: data.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c data.cpp -o data.o

operators.o: operators.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c operators.cpp -o operators.o

linalg.o: linalg.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c linalg.cpp -o linalg.o

main: $(OBJ) main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) *.o main.cpp -o main

clean:
	rm -f main
	rm -f *.o
	rm -f *.i

