CC = gcc -g -O3 -fopenmp
NVCC = nvcc -g -O3 -arch=sm_30
F90 = gfortran -g -O3 -fopenmp

all: sincos

sincos: sincos_host.o sincos_kernel.o cudafor.o cudaLaunch.o main.o
	$(F90) -rdynamic $^ -o $@ -L$(shell dirname $(shell which nvcc))/../lib64 -lcudart -lcuda -ldl

cudaLaunch.o: cudaLaunch.cu
	nvcc -c $< -o $@

sincos_host.o: sincos_host.f90 cudafor.o
	$(F90) -c $< -o $@

cudafor.o: cudafor.f90
	$(F90) -c $< -o $@

sincos_kernel.o: sincos_kernel.cu
	$(NVCC) -c $< -o $@

main.o: main.c
	$(CC) -c $< -o $@

clean:
	rm -rf *.o *.mod sincos

