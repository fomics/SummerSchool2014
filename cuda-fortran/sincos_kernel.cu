// CUDA kernel in C
extern "C" __global__ void sincos_kernel(int nx, int ny, int nz, float* x, float* y, float* xy)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
	
	if ((i >= nx) || (j >= ny) || (k >= nz)) return;

	int index = i + j * nx + k * nx * ny;

	xy[index] = sin(x[index]) + cos(y[index]);
}

