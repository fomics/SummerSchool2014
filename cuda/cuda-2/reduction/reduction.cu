#include <assert.h>
#include <cstdio>

inline __device__ int pow2roundup(int x)
{
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;
    return x;
}

__global__ void reduction_kernel(int x_length, const double* x, int y_length, const double* y, double* result)
{
	extern __shared__ double shared[];

	int half = max(x_length, y_length) / 2 + max(x_length, y_length) % 2;

	// Each block hanldes (2 * blockDim.x) elements of reduction.
	int i = 2 * blockIdx.x * blockDim.x + threadIdx.x;

	// Load product of first 2 pairs into shared memory:
	// idx-th and (idx + blockDim.x)-th.
	shared[threadIdx.x] = 0;
	if (i < x_length)
		shared[threadIdx.x] += x[i];
	if (i < y_length)
		shared[threadIdx.x] += y[i];
	if (i + blockDim.x < x_length)
		shared[threadIdx.x] += x[i + blockDim.x];
	if (i + blockDim.x < y_length)
		shared[threadIdx.x] += y[i + blockDim.x];

	__syncthreads();

	// Reduce pairs in shared memory.
	int start = pow2roundup((blockDim.x + blockDim.x % 2) >> 1);
	//if ((blockDim.x & (blockDim.x - 1)) == 0) start = blockDim.x / 2;
	for (int s = start; s > 32; s >>= 1)
	{
		if ((threadIdx.x < s) && (threadIdx.x + s < half))
			shared[threadIdx.x] += shared[threadIdx.x + s];

		__syncthreads();
	}

	// Unroll last 32 iterations of loop (64 elements).
	// There is no need for synchronizations, since all accesses
	// are within single warp.
	if (threadIdx.x < 32)
	{
		volatile double* vshared = shared;
		if (threadIdx.x + 32 < half) vshared[threadIdx.x] += vshared[threadIdx.x + 32];
		if (threadIdx.x + 16 < half) vshared[threadIdx.x] += vshared[threadIdx.x + 16];
		if (threadIdx.x +  8 < half) vshared[threadIdx.x] += vshared[threadIdx.x +  8];
		if (threadIdx.x +  4 < half) vshared[threadIdx.x] += vshared[threadIdx.x +  4];
		if (threadIdx.x +  2 < half) vshared[threadIdx.x] += vshared[threadIdx.x +  2];
		if (threadIdx.x +  1 < half) vshared[threadIdx.x] += vshared[threadIdx.x +  1];
	}

	// The first thread writes the result.
	if (threadIdx.x == 0)
		result[blockIdx.x] = shared[0];
}

#include <assert.h>
#include <cstdio>
#include <cstdlib>

#define CUDA_ERR_CHECK(x)                                   \
    do { cudaError_t err = x;                               \
        if (err != cudaSuccess) {                           \
        printf("CUDA error %d \"%s\" at %s:%d\n",           \
        (int)err, cudaGetErrorString(err),                  \
        __FILE__, __LINE__); assert(false);                 \
    }} while (0);

#define CUDA_LAUNCH_ERR_CHECK(...)                          \
    do { __VA_ARGS__; cudaError_t err = cudaGetLastError(); \
        if (err != cudaSuccess) {                           \
        printf("CUDA error %d \"%s\" at %s:%d\n",           \
        (int)err, cudaGetErrorString(err),                  \
        __FILE__, __LINE__); assert(false);                 \
    }} while (0);

#include <thrust/extrema.h>

namespace cpu
{
	cudaDeviceProp props;
}

namespace gpu
{
	// We redefine dim3 under namespace, because the default one has
	// constructors, which is not allowed for types device variables
	// (dim3 is used as device vars type below to keep kernel compute
	// grid configuration).
	struct dim3
	{
		unsigned int x, y, z;
		
		__host__ __device__ operator ::dim3()
		{
			return ::dim3(x, y, z);
		}
	};

	struct block_size_to_dynamic_smem_size : public thrust::unary_function<size_t, size_t>
	{
		__host__ __device__
		float operator()(size_t szblock) { return szblock * sizeof(double); }
	};
	
	// Use Thrust occupancy calculator to determine the best size of block.
	template<typename T>
	inline size_t get_optimal_szblock(T kernel)
	{
		using namespace gpu;
		using namespace thrust::system::cuda::detail;

		struct function_attributes_t attrs;
		{
			cudaFuncAttributes funcAttrs;
			CUDA_ERR_CHECK(cudaFuncGetAttributes(&funcAttrs, kernel));
			attrs.constSizeBytes = funcAttrs.constSizeBytes;
			attrs.localSizeBytes = funcAttrs.localSizeBytes;
			attrs.maxThreadsPerBlock = funcAttrs.maxThreadsPerBlock;
			attrs.numRegs = funcAttrs.numRegs;
			attrs.sharedSizeBytes = funcAttrs.sharedSizeBytes;
		}
		struct device_properties_t props;
		{
			props.major = cpu::props.major;
			memcpy(&props.maxGridSize, &cpu::props.maxGridSize, sizeof(int) * 3);
			props.maxThreadsPerBlock = cpu::props.maxThreadsPerBlock;
			props.maxThreadsPerMultiProcessor = cpu::props.maxThreadsPerMultiProcessor;
			props.minor = cpu::props.minor;
			props.multiProcessorCount = cpu::props.multiProcessorCount;
			props.regsPerBlock = cpu::props.regsPerBlock;
			props.sharedMemPerBlock = cpu::props.sharedMemPerBlock;
			props.warpSize = cpu::props.warpSize;
		}
		return block_size_with_maximum_potential_occupancy(attrs, props, block_size_to_dynamic_smem_size());
	}

	template<typename T>
	inline void get_optimal_grid_block_config(T kernel,
		int nx, int ny, dim3* grid, dim3* blocks)
	{
		size_t szblock = get_optimal_szblock(kernel);

		grid->x = 1; grid->y = 1; grid->z = 1;
		blocks->x = 1; blocks->y = 1; blocks->z = 1;

		if (szblock > nx)
		{
			blocks->x = nx;
			blocks->y = min(ny, (int)szblock / blocks->x);
			grid->y = ny / blocks->y;
			if (ny % blocks->y) grid->y++;
		}
		else
		{
			blocks->x = szblock;
			grid->x = nx / blocks->x;
			if (nx % blocks->x) grid->x++;
			grid->y = ny;
		}
	}

	#define determine_optimal_grid_block_config(kernel_name, nx, ny) \
	{ \
		{ \
			gpu::config_t c; \
			gpu::get_optimal_grid_block_config(kernel_name, nx, ny, &c.grid, &c.block); \
			memcpy(&gpu::config, &c, sizeof(gpu::config_t)); \
		} \
	}

	typedef struct __attribute__((packed)) { dim3 grid, block; } config_t;

	config_t config;
}

__global__ void check()
{
	assert(warpSize == 32);
}

double reduction(int length, double* in1, double* in2)
{
	using namespace gpu;

	check<<<1, 1>>>();

	determine_optimal_grid_block_config(reduction_kernel, length / 2 + length % 2, 1);
	
	double* buffer = NULL;
	CUDA_ERR_CHECK(cudaMalloc(&buffer, sizeof(double) * (config.grid.x + config.grid.x % 2)));

	double* x_dev = in1;
	double* y_dev = in2;
	for (int szbuffer = config.grid.x, x_length = length, y_length = length; ; szbuffer = config.grid.x)
	{
		CUDA_LAUNCH_ERR_CHECK(reduction_kernel<<<config.grid, config.block, config.block.x * sizeof(double)>>>(
			x_length, x_dev, y_length, y_dev, buffer));

		if (szbuffer == 1) break;

		x_length = szbuffer / 2 + szbuffer % 2;
		determine_optimal_grid_block_config(reduction_kernel, x_length / 2 + x_length % 2, 1);
		y_length = szbuffer / 2;

		x_dev = buffer;
		y_dev = buffer + x_length;
	}

	double result;
	CUDA_ERR_CHECK(cudaMemcpy(&result, &buffer[0], sizeof(double), cudaMemcpyDeviceToHost));

	CUDA_ERR_CHECK(cudaFree(buffer));
	
	return result;
}

int main(int argc, char* argv[])
{
	using namespace gpu;

	if (argc != 2)
	{
		printf("Usage: %s <length>\n", argv[0]);
		return 1;
	}
	
	int length = atoi(argv[1]);

	double *x = (double*)malloc(sizeof(double) * length);
	double *y = (double*)malloc(sizeof(double) * length);
#if 0
	for (int i = 0; i < length; i++)
	{
		x[i] = 2 * i + 1; /*rand();*/ printf("%d ", 2 * i + 1);
		y[i] = 2 * i + 2; /*rand();*/ printf("%d ", 2 * i + 2);
	}
	printf("\n");
#else
	for (int i = 0; i < length; i++)
	{
		x[i] = rand() % 2;
		y[i] = rand() % 2;
	}
#endif

	CUDA_ERR_CHECK(cudaGetDeviceProperties(&cpu::props, 0));

	if (length == 1) return 0; // TODO

	CUDA_ERR_CHECK(cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024 * 1024 * 32));

	double *x_dev; CUDA_ERR_CHECK(cudaMalloc(&x_dev, sizeof(double) * length));
	double *y_dev; CUDA_ERR_CHECK(cudaMalloc(&y_dev, sizeof(double) * length));

	CUDA_ERR_CHECK(cudaMemcpy(x_dev, x, sizeof(double) * length, cudaMemcpyHostToDevice));
	CUDA_ERR_CHECK(cudaMemcpy(y_dev, y, sizeof(double) * length, cudaMemcpyHostToDevice));

	double result_gpu = reduction(length, x_dev, y_dev);

	CUDA_ERR_CHECK(cudaFree(x_dev));
	CUDA_ERR_CHECK(cudaFree(y_dev));
	
	length = atoi(argv[1]);
	double result_cpu = 0;
	for (int i = 0; i < length; i++)
		result_cpu += x[i] + y[i];

	printf("GPU result = %f\n", result_gpu);
	printf("CPU result = %f\n", result_cpu);
	
	if (fabs(result_gpu - result_cpu) > 1e-6)
		printf("FAIL\n");
	else
		printf("OK!\n");
	
    return 0;
}

