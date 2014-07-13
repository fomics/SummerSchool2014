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

struct SoA
{
    uint8_t *r;
    uint8_t *g;
    uint8_t *b;
};

__global__ void kernel_SOA(uint8_t* r, uint8_t* g, uint8_t* b, int N)
{
	int idxx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idxx < N)
	{
		uint8_t b1, b2, b3, res;
		b1 = r[idxx];
		b2 = g[idxx];
		b3 = b[idxx];
		res = 0.2126f * b1 + 0.7152f * b2 + 0.0722f * b3;
		r[idxx] = res;
		g[idxx] = res;
		b[idxx] = res;
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
	size_t size;
	dim3 block;
	block.x = BLOCK_SIZE;
	dim3 grid(((N - 1) / block.x + 1), 1, 1);

	size = N * sizeof(uint8_t);

	SoA* hSoA;
	hSoA = new SoA();
	hSoA->r = (uint8_t*)malloc(size);
	hSoA->g = (uint8_t*)malloc(size);
	hSoA->b = (uint8_t*)malloc(size);
	
	uint8_t *r, *g, *b;
	CUDA_ERR_CHECK(cudaMalloc((void**)&r,size));
	CUDA_ERR_CHECK(cudaMalloc((void**)&g,size));
	CUDA_ERR_CHECK(cudaMalloc((void**)&b,size));
	
	CUDA_ERR_CHECK(cudaMemcpy(r,hSoA->r,size, cudaMemcpyHostToDevice));
	CUDA_ERR_CHECK(cudaMemcpy(g,hSoA->g,size, cudaMemcpyHostToDevice));
	CUDA_ERR_CHECK(cudaMemcpy(b,hSoA->b,size, cudaMemcpyHostToDevice));
	
	CUDA_ERR_CHECK(cudaEventRecord(start));
	
	kernel_SOA<<<grid,block>>>(r, g, b, N);
	CUDA_ERR_CHECK(cudaGetLastError());
	
	CUDA_ERR_CHECK(cudaEventRecord(stop));
	CUDA_ERR_CHECK(cudaDeviceSynchronize());
	CUDA_ERR_CHECK(cudaEventElapsedTime(&time, start, stop));
	
	CUDA_ERR_CHECK(cudaMemcpy(hSoA->r,r,size, cudaMemcpyDeviceToHost));
	CUDA_ERR_CHECK(cudaMemcpy(hSoA->g,g,size, cudaMemcpyDeviceToHost));
	CUDA_ERR_CHECK(cudaMemcpy(hSoA->b,b,size, cudaMemcpyDeviceToHost));
	
	cout << endl << time << endl;
	
	CUDA_ERR_CHECK(cudaFree(r));
	CUDA_ERR_CHECK(cudaFree(g));
	CUDA_ERR_CHECK(cudaFree(b));

	CUDA_ERR_CHECK(cudaEventDestroy(start));
	CUDA_ERR_CHECK(cudaEventDestroy(stop));

	free(hSoA);
	
	return 0;
}

