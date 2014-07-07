#include <cstdio>
#include <cstdlib>
#include <time.h>

#ifndef CPU
#include <cuda_runtime.h>
#include <cublas_v2.h>
#define CUDA_SAFE_CALL(x) \
	do { cudaError_t err = x; if (err != cudaSuccess) { \
		fprintf (stderr, "Error \"%s\" at %s:%d \n", cudaGetErrorString(err), \
		__FILE__, __LINE__); exit(EXIT_FAILURE); \
	}} while (0);
#define CUBLAS_SAFE_CALL(x) \
	do { cublasStatus_t err = x; if (err != CUBLAS_STATUS_SUCCESS) { \
		fprintf (stderr, "Error at %s:%d \n", \
		__FILE__, __LINE__); exit(EXIT_FAILURE); \
	}} while (0);
#else
#include <mkl.h>
#endif

#if defined(HAVE_SINGLE) && defined(HAVE_DOUBLE)
#error "Macros HAVE_SINGLE and HAVE_DOUBLE cannot be defined simultaneously"
#endif
#ifdef HAVE_SINGLE
#define real float
#define cpu_gemm cblas_sgemm
#define gpu_gemm cublasSgemm
#elif HAVE_DOUBLE
#define real double
#define cpu_gemm cblas_dgemm
#define gpu_gemm cublasDgemm
#else
#error "Macro HAVE_SINGLE or HAVE_DOUBLE must be defined"
#endif

// Get the timer float.
static void get_time(double* ret)
{
	volatile struct timespec val;
	clock_gettime(CLOCK_REALTIME, (struct timespec*)&val);
	*ret = (double)0.000000001 * val.tv_nsec + val.tv_sec;
}

void usage(const char* filename)
{
	printf("Dense n x n matrix multiply: C := A * B.\n");
	printf("Usage: %s <n>\n", filename);
}

// Memory alignment, for vectorization on MIC.
// 4096 should be best for memory transfers over PCI-E.
#define MEMALIGN 4096

int main(int argc, char* argv[])
{
	const int printable_n = 128;

	if (argc != 2)
	{
		usage(argv[0]);
		return 0;
	}

	int n = atoi(argv[1]);
	if (n <= 0)
	{
		usage(argv[0]);
		return 0;
	}

	real* A; posix_memalign((void**)&A, MEMALIGN, n * n * sizeof(real));
	real* B; posix_memalign((void**)&B, MEMALIGN, n * n * sizeof(real));
	real* C; posix_memalign((void**)&C, MEMALIGN, n * n * sizeof(real));

	// Generate the random dense matrices.
	for (int i = 0; i < n * n; i++)
	{
		A[i] = drand48();
		B[i] = drand48();
	}

	// Print out the input data if n is small.
	if (n <= printable_n)
	{
		printf("Input data:\n");
		for (int i = 0; i < n; i++)
			printf("(%f, %f)\n", A[i], B[i]);
		printf("\n");
	}

	// Strip initialization delay.
	double init_start, init_finish;
	get_time(&init_start);
#ifndef CPU
	int count = 0;
	cudaGetDeviceCount(&count);
#endif
	get_time(&init_finish);

	double load_start, load_finish;
	get_time(&load_start);
#ifndef CPU
	real *A_dev = NULL, *B_dev = NULL, *C_dev = NULL;
	CUDA_SAFE_CALL(cudaMalloc(&A_dev, n * n * sizeof(real)));
	CUDA_SAFE_CALL(cudaMalloc(&B_dev, n * n * sizeof(real)));
	CUDA_SAFE_CALL(cudaMalloc(&C_dev, n * n * sizeof(real)));

	cublasHandle_t handle;	
	CUBLAS_SAFE_CALL(cublasCreate(&handle));

	CUBLAS_SAFE_CALL(cublasSetMatrix(n, n, sizeof(real), A, n, A_dev, n));
	CUBLAS_SAFE_CALL(cublasSetMatrix(n, n, sizeof(real), B, n, B_dev, n));
	CUBLAS_SAFE_CALL(cublasSetMatrix(n, n, sizeof(real), C, n, C_dev, n));
#endif
	get_time(&load_finish);

	double compute_start, compute_finish;
	get_time(&compute_start);
#ifndef CPU
	real alpha = 1.0, beta = 0.0;
	CUBLAS_SAFE_CALL(gpu_gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
		n, n, n, &alpha, A_dev, n, B_dev, n, &beta, C_dev, n));
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
#else
	real alpha = 1.0, beta = 0.0;
	cpu_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		n, n, n, alpha, A, n, B, n, beta, C, n);
#endif
	get_time(&compute_finish);

	double save_start, save_finish;
	get_time(&save_start);
#ifndef CPU
	CUBLAS_SAFE_CALL(cublasGetMatrix(n, n, sizeof(real), C_dev, n, C, n));

	CUDA_SAFE_CALL(cudaFree(A_dev));
	CUDA_SAFE_CALL(cudaFree(B_dev));
	CUDA_SAFE_CALL(cudaFree(C_dev));
	cublasDestroy(handle);
#endif
	get_time(&save_finish);

	printf("Init time = %f sec\n", init_finish - init_start);
	printf("Load time = %f sec\n", load_finish - load_start);
	printf("Compute time = %f sec ~ %f GFLOPS\n", compute_finish - compute_start,
		2 * 1.0e-9 * n * n * n / (compute_finish - compute_start));
	printf("Save time = %f sec\n", save_finish - save_start);

	// Print out the output data if n is small.
	if (n <= printable_n)
	{
		printf("Output data:\n");
		for (int i = 0; i < n; i++)
			printf("%f\n", C[i]);
		printf("\n");
	}
#if 1
	// Show last 10 pairs.
	printf("Last 10 pairs:\n");
	for (int i = n - 11; i < n; i++)
	{
		printf("%f\n", C[i]);
	}
	printf("\n");
#endif

	free(A);
	free(B);
	free(C);

	return 0;
}

