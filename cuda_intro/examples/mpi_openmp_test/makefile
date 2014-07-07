all: mpi_openmp_test

mpi_openmp_test: mpi_openmp_test.c
	mpicc -O3 -fopenmp $< -o $@

clean:
	rm -rf mpi_openmp_test

test: mpi_openmp_test
	OMP_NUM_THREADS=2 mpirun -np 2 ./mpi_openmp_test

