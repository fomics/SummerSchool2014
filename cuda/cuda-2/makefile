ARCH = 30
NVCC = nvcc
CXXFLAGS = -mtune=native -ffast-math -fopenmp
NVFLAGS = -Icuda_launch_config -g -O3 -rdc=true -lineinfo -arch=sm_$(ARCH) -use_fast_math $(addprefix -Xcompiler=,$(CXXFLAGS)) -Xptxas=-v -lcudadevrt

HEADERS = stats.h   data.h   operators.h   linalg.h
OBJ     = stats.o   data.o   main.o

.SUFFIXES: .cu

all: main

stats.o: stats.cu stats.h
	$(NVCC) $(NVFLAGS) $(IFLAGS) -c $< -o $@

data.o: data.cu data.h
	$(NVCC) $(NVFLAGS) $(IFLAGS) -c $< -o $@

main.o: main.cu $(HEADERS)
	$(NVCC) $(NVFLAGS) $(IFLAGS) -c $< -o $@

main: $(OBJ)
	$(NVCC) $(NVFLAGS) $^ -o $@

clean:
	rm -f main
	rm -f *.o
	rm -f *.i

test: main
	./main 128 128 100 0.01 > result && diff -u result result.ref

