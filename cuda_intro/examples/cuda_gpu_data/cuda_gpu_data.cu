#include <stdio.h>

#define CUDA_ERR_CHECK(x) \
	do { cudaError_t err = x; if (err != cudaSuccess) { \
		fprintf (stderr, "Error \"%s\" at %s:%d \n", \
		 cudaGetErrorString(err), \
		__FILE__, __LINE__); exit(-1); \
	}} while (0);

__global__ void gpu_kernel(int* gpu_data)
{
	int4 coords = { blockIdx.x, gridDim.x, threadIdx.x, blockDim.x };
	((int4*)gpu_data)[blockIdx.x * blockDim.x + threadIdx.x] = coords;
}

int main()
{
	int nthreads = 4, nblocks = 4;
	int n = nthreads * nblocks * 4;
	int* host_data = (int*)malloc(n * sizeof(int));
	int* gpu_data = NULL;
	CUDA_ERR_CHECK( cudaMalloc(&gpu_data, n * sizeof(int)) );

	gpu_kernel<<<nblocks, nthreads>>>(gpu_data);

	CUDA_ERR_CHECK( cudaMemcpy(host_data, gpu_data, n * sizeof(int),
		cudaMemcpyDeviceToHost) );
	CUDA_ERR_CHECK( cudaFree(gpu_data) );
	for (int i = 0; i < n; i += 4)
		printf("Hello form b #%d of %d, t #%d of %d!\n",
			host_data[i], host_data[i+1], host_data[i+2], host_data[i+3]);

	free(host_data);
	
	return 0;
}

