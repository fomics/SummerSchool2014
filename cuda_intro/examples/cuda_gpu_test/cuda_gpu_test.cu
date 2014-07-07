#include <stdio.h>

__global__ void gpu_kernel()
{
	int block_idx, grid_dim;
	block_idx = blockIdx.x;
	grid_dim = gridDim.x;
	   
	printf("Hello form b #%d of %d, t #%d of %d!\n",
		block_idx, grid_dim, threadIdx.x, blockDim.x);			
}

int main()
{
	gpu_kernel<<<2, 2>>>();
	cudaDeviceSynchronize();
	
	return 0;
}

