CXX= cc
CXXFLAGS=-O3
CXXFLAGS=-O3 -h acc

SOURCES = stats.c data.c operators.c linalg.c
HEADERS = stats.h   data.h   operators.h   linalg.h
OBJ     = stats.o   data.o   operators.o   linalg.o

.SUFFIXES: .c

all: main

stats.o: stats.c stats.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c stats.c -o stats.o

data.o: data.c data.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c data.c -o data.o

operators.o: operators.c operators.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c operators.c -o operators.o

linalg.o: linalg.c linalg.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c linalg.c -o linalg.o

main: $(OBJ) main.c $(HEADERS)
	$(CXX) $(CXXFLAGS) *.o main.c -o main

clean:
	rm -f main
	rm -f *.o
	rm -f *.i

