#include <mpi.h>
#include <omp.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	
	int block_idx, grid_dim;
	MPI_Comm_rank(MPI_COMM_WORLD, &block_idx);
	MPI_Comm_size(MPI_COMM_WORLD, &grid_dim);
	
	#pragma omp parallel
	{
		printf("Hello form b #%d of %d, t #%d of %d!\n",
			block_idx, grid_dim,
			omp_get_thread_num(), omp_get_num_threads());
	}
	MPI_Finalize();
	return 0;
}
