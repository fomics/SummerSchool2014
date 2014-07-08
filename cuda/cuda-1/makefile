NVCC = nvcc -arch=sm_35
CXX = gcc
CXXFLAGS = -O3 -g -fopenmp

HEADERS = stats.h   data.h   operators.h   linalg.h
OBJ     = stats.o   data.o   operators.o   linalg.o   main.o

all: main

stats.o: stats.c stats.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c stats.c -o stats.o

data.o: data.c data.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c data.c -o data.o

operators.o: operators.c operators.h
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c operators.c -o operators.o

linalg.o: linalg.cu linalg.h
	$(NVCC) $(addprefix -Xcompiler , $(CXXFLAGS)) -c linalg.cu -o linalg.o

main.o: main.c $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

main: $(OBJ) $(HEADERS)
	$(NVCC) -Xcompiler -fopenmp *.o -o main -lm

clean:
	rm -f main
	rm -f *.o
	rm -f *.i

