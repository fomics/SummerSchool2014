#include <iostream>
#include <stdio.h>
#include <stdint.h>

#define CUDA_ERR_CHECK(x) \
        do { cudaError_t err = x; if (err != cudaSuccess) { \
                fprintf (stderr, "Error \"%s\" at %s:%d \n", \
                 cudaGetErrorString(err), \
                __FILE__, __LINE__); exit(-1); \
        }} while (0);

#define BLOCK_SIZE 256

struct AoS
{
    uint8_t r, g, b;
};

__global__ void kernel_AOS(AoS* A, int N)
{
	int idxx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idxx < N)
	{
		uint8_t b1, b2, b3, res;
		b1 = A[idxx].r;
		b2 = A[idxx].g;
		b3 = A[idxx].b;
		res = 0.2126f * b1 + 0.7152f * b2 + 0.0722f * b3;
		A[idxx].r = res;
		A[idxx].g = res;
		A[idxx].b = res;
	}
}

int main()
{
	using namespace std;

	int N;
	cout << "input N" << endl;
	cin >> N;
	
	cudaEvent_t start, stop;
	CUDA_ERR_CHECK(cudaEventCreate(&start));
	CUDA_ERR_CHECK(cudaEventCreate(&stop));
	
	float time;
	size_t size = N * sizeof(AoS);
	
	AoS *hAoS, *dAoS;
	hAoS = (AoS*)malloc(size);
	CUDA_ERR_CHECK(cudaMalloc((void**)&dAoS, size));
	
	dim3 block;
	block.x = BLOCK_SIZE;
	dim3 grid(((N - 1) / block.x + 1), 1, 1);
	
	CUDA_ERR_CHECK(cudaMemcpy(dAoS, hAoS, size, cudaMemcpyHostToDevice));
	CUDA_ERR_CHECK(cudaEventRecord(start));
	
	kernel_AOS<<<grid, block>>>(dAoS, N);
	CUDA_ERR_CHECK(cudaGetLastError());
	
	CUDA_ERR_CHECK(cudaEventRecord(stop));
	CUDA_ERR_CHECK(cudaDeviceSynchronize());
	CUDA_ERR_CHECK(cudaEventElapsedTime(&time, start, stop));
	CUDA_ERR_CHECK(cudaMemcpy(hAoS,dAoS,size, cudaMemcpyDeviceToHost));
	
	cout << endl << time << endl;
	CUDA_ERR_CHECK(cudaFree(dAoS));
	free(hAoS);

	CUDA_ERR_CHECK(cudaEventDestroy(start));
	CUDA_ERR_CHECK(cudaEventDestroy(stop));
	
	return 0;
}

