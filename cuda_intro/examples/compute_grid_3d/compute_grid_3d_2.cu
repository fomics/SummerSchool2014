#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define CUDA_ERR_CHECK(x) \
	do { cudaError_t err = x; if (err != cudaSuccess) { \
		fprintf (stderr, "Error \"%s\" at %s:%d \n", \
		cudaGetErrorString(err), \
		__FILE__, __LINE__); exit(-1); \
	}} while (0);

#define roundup(n, width) (((n) + (width) - 1) & ~unsigned((width) - 1))
#define rounddown(n, width) (n / width * width)

__global__ void gpu_kernel(int ni, int nj, int nk, int* data, int i_inc, int j_inc, int k_inc)
{
	for (int k = blockIdx.z * blockDim.z + threadIdx.z; k < nk; k += k_inc)
		for (int j = blockIdx.y * blockDim.y + threadIdx.y; j < nj; j += j_inc)
			for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < ni; i += i_inc)
			{
				int idx = i + ni * j + nj * ni * k;
				assert(data[idx] == 0);
				data[idx] = idx;
			}
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("Usage: %s <ni> <nj> <nk>\n", argv[0]);
		return 0;
	}

	int ni = atoi(argv[1]);
	int nj = atoi(argv[2]);
	int nk = atoi(argv[3]);

	int* data = (int*)malloc(ni * nj * nk * sizeof(int));
	for (int k = 0; k < nk; k++)
		for (int j = 0; j < nj; j++)
			for (int i = 0; i < ni; i++)
			{
				int idx = i + ni * j + nj * ni * k;
				data[idx] = idx;
			}

	int* gpu_data;
	CUDA_ERR_CHECK( cudaMalloc(&gpu_data, ni * nj * nk * sizeof(int)) );
	CUDA_ERR_CHECK( cudaMemset(gpu_data, 0, ni * nj * nk * sizeof(int)) );

	struct cudaDeviceProp props;
	CUDA_ERR_CHECK( cudaGetDeviceProperties(&props, 0) );
	dim3 max_grid;
	max_grid.x = props.maxGridSize[0];
	max_grid.y = props.maxGridSize[1];
	max_grid.z = props.maxGridSize[2];
	
	gpu_kernel<<<dim3(min(max(1, roundup(ni, 16) / 16), rounddown(max_grid.x, 16)),
	                  min(max(1, roundup(nj,  8) /  8), rounddown(max_grid.y,  8)),
	                  min(max(1, roundup(nk,  8) /  8), rounddown(max_grid.z,  8))),
	             dim3(16, 8, 8)>>>(ni, nj, nk, gpu_data,
	                  min(max(1, roundup(ni, 16) / 16), rounddown(max_grid.x, 16)) * 16,
	                  min(max(1, roundup(nj,  8) /  8), rounddown(max_grid.y,  8)) * 8,
	                  min(max(1, roundup(nk,  8) /  8), rounddown(max_grid.z,  8)) * 8);

	CUDA_ERR_CHECK( cudaGetLastError() );
	
	int* host_data = (int*)malloc(ni * nj * nk * sizeof(int));
	CUDA_ERR_CHECK( cudaMemcpy(host_data, gpu_data, ni * nj * nk * sizeof(int),
		cudaMemcpyDeviceToHost) );
	CUDA_ERR_CHECK( cudaFree(gpu_data) );
	
	for (int k = 0; k < nk; k++)
		for (int j = 0; j < nj; j++)
			for (int i = 0; i < ni; i++)
			{
				int idx = i + ni * j + nj * ni * k;
				if (data[idx] != host_data[idx])
				{
					fprintf(stderr, "Values mismatch at (i, j, k) = (%d, %d, %d): %d != %d\n",
						i, j, k, data[idx], host_data[idx]);
					exit(1);
				}
			}

	printf("Done!\n");

	free(host_data);
	free(data);

	return 0;
}
