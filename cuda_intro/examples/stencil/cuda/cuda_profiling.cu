//===----------------------------------------------------------------------===//
//
//     KernelGen -- A prototype of LLVM-based auto-parallelizing Fortran/C
//        compiler for NVIDIA GPUs, targeting numerical modeling code.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include <dlfcn.h>
#include <stdarg.h>
#include <stdio.h>

using namespace std;

#include "cuda_profiling.h"
#include "timing.h"

static char* wrapper_funcname = 0;
__attribute__((unused)) static long wrapper_lineno = 0;

int kernel_enable_regcount(char* funcname, long lineno)
{
	wrapper_funcname = funcname;
	wrapper_lineno = lineno;
	return 0;
}

int kernel_disable_regcount()
{
	wrapper_funcname = 0;
	return 0;
}

int kernel_configure_gird(int ndims, int nx, int ny, int ns,
	kernel_config_t* config)
{
	config->szshmem = 0;
	dim3* blockDim = &config->blockDim;
	dim3* gridDim = &config->gridDim;
	dim3* strideDim = &config->strideDim;

	// For 1D case block dimensions are fixed to the best values we've
	// found earlier during exhaustive search within KernelGen.
	switch (ndims)
	{
	case 1 :
		blockDim->x = 128;
		blockDim->y = 1;
		blockDim->z = 1;
		break;
	case 2 :
		blockDim->x = 32;
		blockDim->y = 16;
		blockDim->z = 1;
		break;
	default :
		fprintf(stderr, "Unsupported number of dimensions: %d\n", ndims);
		return -1;
	}
	
	strideDim->x = nx;
	strideDim->y = ny;
	strideDim->z = ns;
	
	// If maximum permitted GPU grid dimensions are large
	// enough to hold the original problem dimensions,
	// we just divide the problem size by block sizes.
	// Otherwise, the only solution is to assign multiple
	// grid points to each block thread. In order to comply
	// with coalescing requirements, we make a gap (stride)
	// between grid points indexes assigned to the same thread.
#if defined(__CUDA_VECTORIZE2__)
	gridDim->x = nx / (2 * blockDim->x);
#elif defined(__CUDA_VECTORIZE4__)
	gridDim->x = nx / (4 * blockDim->x);
#else
	gridDim->x = nx / blockDim->x;
#endif
	gridDim->y = ny / blockDim->y;
	gridDim->z = ns / blockDim->z;
#if defined(__CUDA_VECTORIZE2__)
	if (nx % (2 * blockDim->x)) gridDim->x++;
#elif defined(__CUDA_VECTORIZE4__)
	if (nx % (4 * blockDim->x)) gridDim->x++;
#else
	if (nx % blockDim->x) gridDim->x++;
#endif
	if (ny % blockDim->y) gridDim->y++;
	if (ns % blockDim->z) gridDim->z++;
	struct cudaDeviceProp props;
	CUDA_SAFE_CALL(cudaGetDeviceProperties(&props, 0));
	if (props.maxGridSize[0] * blockDim->x < nx)
	{
		gridDim->x = props.maxGridSize[0];
		strideDim->x = props.maxGridSize[0] * blockDim->x;
	}
	if (props.maxGridSize[1] * blockDim->y < ny)
	{
		gridDim->y = props.maxGridSize[1];
		strideDim->y = props.maxGridSize[1] * blockDim->y;
	}
	if (props.maxGridSize[2] * blockDim->z < ns)
	{
		gridDim->z = props.maxGridSize[2];
		strideDim->z = props.maxGridSize[2] * blockDim->z;
	}
	
	return 0;
}

// TODO We may need to combine this call with configure_grid^^ as maximum shared
// memory size may influence the grid config.
int kernel_configure_shmem(kernel_config_t* config,
	int narrays, ...)
{
	// We use shared memory as a scratch space for threads block.
	ptrdiff_t* shmem_arrays = (ptrdiff_t*)malloc(sizeof(ptrdiff_t) * narrays);
	va_list list;
	va_start(list, narrays);
	size_t szshmem = 0;
	for (int i = 0; i < narrays; i++)
	{
		size_t szelem = va_arg(list, size_t);
		size_t num_block_elems = va_arg(list, size_t);
		size_t offset = va_arg(list, size_t);
#if defined(__CUDA_VECTORIZE2__)
		shmem_arrays[i] = szshmem + offset * szelem * 2;
		szshmem += num_block_elems * szelem * 2;
#elif defined(__CUDA_VECTORIZE4__)
		shmem_arrays[i] = szshmem + offset * szelem * 4;
		szshmem += num_block_elems * szelem * 4;
#else
		shmem_arrays[i] = szshmem + offset * szelem;
		szshmem += num_block_elems * szelem;
#endif
	}
	va_end(list);

	// Make sure shared memory size does not exceed the maximum amount
	// allowed for the threads block.
	struct cudaDeviceProp prop;
	CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, 0));
	size_t szshmem_max = prop.sharedMemPerBlock;
	if (szshmem > szshmem_max)
	{
		fprintf(stderr, "Required amount of shared memory %zu exceeds the "
			"maximum allowed for one threads block (%zu)\n",
			szshmem, szshmem_max);
		return -1;
	}

	config->szshmem = szshmem;
	CUDA_SAFE_CALL(cudaMalloc(&config->shmem_arrays,
		sizeof(ptrdiff_t) * narrays));
	CUDA_SAFE_CALL(cudaMemcpy(config->shmem_arrays,
		shmem_arrays, sizeof(ptrdiff_t) * narrays, cudaMemcpyHostToDevice));
	free(shmem_arrays);
	
	return 0;
}

void kernel_config_dispose(kernel_config_t* config)
{
	CUDA_SAFE_CALL(cudaFree(config->shmem_arrays));
}

extern "C" cudaError_t __real_cudaLaunch(const void *func);

extern "C" cudaError_t __wrap_cudaLaunch(const void *func)
{
	if (!wrapper_funcname)
		return __real_cudaLaunch(func);

	// Find out the kernel name.
	Dl_info info;
	if (dladdr(func, &info) == 0)
	{
		fprintf(stderr, "Error in dladdr(%p, %p): %s\n", func, &info, dlerror());
		exit(-1);
	}
	if (info.dli_saddr != func)
	{
		fprintf(stderr, "Cannot find kernel name for address %p\n", func);
		exit(-1);
	}
	const char* kernel_name = info.dli_sname;

	// Get the kernel register count.
	struct cudaFuncAttributes attrs;
	CUDA_SAFE_CALL(cudaFuncGetAttributes(&attrs, func));
	printf("%s regcount = %d\n", kernel_name, attrs.numRegs);

	// Get the kernel execution time.
	struct timespec start, finish;
	get_time(&start);
	cudaError_t result = __real_cudaLaunch(func);
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
	get_time(&finish);
	double kernel_time = get_time_diff(&start, &finish);
	if (kernel_name)
		printf("%s kernel time = %f\n", kernel_name, kernel_time);
	else
		printf("kernel time = %f\n", kernel_time);
	return result;
}

