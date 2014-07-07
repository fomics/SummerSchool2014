//===----------------------------------------------------------------------===//
//
//     KernelGen -- A prototype of LLVM-based auto-parallelizing Fortran/C
//        compiler for NVIDIA GPUs, targeting numerical modeling code.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#define CUDA_SAFE_CALL(x) \
	do { cudaError_t err = x; if (err != cudaSuccess) { \
		fprintf (stderr, "Error \"%s\" at %s:%d \n", cudaGetErrorString(err), \
		__FILE__, __LINE__); exit(-1); \
	}} while (0);

#ifdef __cplusplus
extern "C" {
#endif

int kernel_enable_regcount(char* funcname, long lineno);

int kernel_disable_regcount();

struct kernel_config_t
{
	dim3 gridDim, blockDim, strideDim;
	
	// The total size of shared memory requried
	// (to be used with the kernel launch call).
	size_t szshmem;
	
	// The array of offsets where the corresponding
	// shared memory scratch buffers of input data
	// arrays shall start.
	ptrdiff_t* shmem_arrays;
};

int kernel_configure_gird(int ndims, int nx, int ny, int ns,
	kernel_config_t* config);

int kernel_configure_shmem(kernel_config_t* config,
	int narrays, ...);

void kernel_config_dispose(kernel_config_t* config);

#ifdef __cplusplus
}
#endif

