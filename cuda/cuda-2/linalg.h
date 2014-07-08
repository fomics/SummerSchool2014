// linear algebra subroutines
// Ben Cumming @ CSCS

#ifndef LINALG_H
#define LINALG_H

#include "check.h"
#include "data.h"
#include "linalg.h"
#include "operators.h"
#include "stats.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace
{
	int cg_initialized = 0;
	double *r, *Ap, *p;
	double *Fx, *Fxold, *v, *xold; // 1d

	// initialize temporary storage fields used by the cg solver
	// I do this here so that the fields are persistent between calls
	// to the CG solver. This is useful if we want to avoid malloc/free calls
	// on the device for the OpenACC implementation (feel free to suggest a better
	// method for doing this)
	inline void cg_init(const int N)
	{
		using namespace gpu;
	
		CUDA_ERR_CHECK(cudaMalloc(&Ap,	sizeof(double) * N + (1 << 7)));
		CUDA_ERR_CHECK(cudaMalloc(&r,	 sizeof(double) * N + (1 << 7)));
		CUDA_ERR_CHECK(cudaMalloc(&p,	 sizeof(double) * N + (1 << 7)));
		CUDA_ERR_CHECK(cudaMalloc(&Fx,	sizeof(double) * N + (1 << 7)));
		CUDA_ERR_CHECK(cudaMalloc(&Fxold, sizeof(double) * N + (1 << 7)));
		CUDA_ERR_CHECK(cudaMalloc(&v,	 sizeof(double) * N + (1 << 7)));
		CUDA_ERR_CHECK(cudaMalloc(&xold,  sizeof(double) * N + (1 << 7)));
	
		cg_initialized = 1;
	}

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
}

////////////////////////////////////////////////////////////////////////////////
//  blas level 1 reductions
////////////////////////////////////////////////////////////////////////////////

namespace gpu
{
	namespace ss_sum_kernel
	{
		// computes the sum of x and y
		// x and y are vectors of lenghts x_length and y_length
		template<short V, typename T>
		__global__ void kernel(
			const int x_length, const double* const __restrict__ x,
			const int y_length, const double* const __restrict__ y, double* __restrict__ result)
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
			for (int s = start; s > warpSize; s >>= 1)
			{
				if ((threadIdx.x < s) && (threadIdx.x + s < half))
					shared[threadIdx.x] += shared[threadIdx.x + s];

				__syncthreads();
			}

			// Unroll last 32 iterations of loop (64 elements).
			// There is no need for synchronizations, since all accesses
			// are within single warp.
			if (threadIdx.x < warpSize)
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

		config_t configs[MAX_CONFIGS];
		double* buffer = NULL;
		size_t szbuffer = 0;
	}
}

// computes the sum of x and y elements
// x and y are vectors of length N
inline double ss_sum(
	const double* const __restrict__ x, const double* const __restrict__ y, const int length)
{
	using namespace gpu;
	using namespace gpu::ss_sum_kernel;

	{
		size_t size = configs[0].grid.x + configs[0].grid.x % 2;
		if (szbuffer != size)
		{
			if (szbuffer)
				CUDA_ERR_CHECK(cudaFree(buffer));
			CUDA_ERR_CHECK(cudaMalloc(&buffer, sizeof(double) * size));
			szbuffer = size;
		}

		CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<
			configs[0].grid, configs[0].block, configs[0].block.x * sizeof(double)>>>(
			length, x, length, y, buffer));
	}

	for (int i = 1, szbuffer = configs[0].grid.x; szbuffer != 1; i++)
	{
		int x_length = szbuffer / 2 + szbuffer % 2;
		int y_length = szbuffer / 2;

		const double* x_dev = buffer;
		const double* y_dev = buffer + x_length;

		CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<
			configs[i].grid, configs[i].block, configs[i].block.x * sizeof(double)>>>(
			x_length, x_dev, y_length, y_dev, buffer));

		szbuffer = configs[i].grid.x;
	}

	// record the number of floating point oporations
	flops_blas1 += length;

	CUDA_ERR_CHECK(cudaDeviceSynchronize());

	double result;
	CUDA_ERR_CHECK(cudaMemcpy(&result, &buffer[0], sizeof(double), cudaMemcpyDeviceToHost));
	return result;
}

// computes the sum of x elements
// x is a vector of length N
inline double ss_sum(const double* const __restrict__ x, const int N)
{
	return ss_sum(x, x + N / 2 + N % 2, N / 2 + N % 2);
}

namespace gpu
{
	namespace ss_dot_kernel
	{
		// computes the inner product of x and y
		// x and y are vectors of length N
		template<short V, typename T>
		__global__ void kernel(const int length,const double* const __restrict__ x,
			const double* const __restrict__ y, double* __restrict__ result)
		{
			extern __shared__ double shared[];

			int half = length / 2 + length % 2;

			// Each block hanldes (2 * blockDim.x) elements of reduction.
			int i = 2 * blockIdx.x * blockDim.x + threadIdx.x;

			// Load product of first 2 pairs into shared memory:
			// idx-th and (idx + blockDim.x)-th.
			shared[threadIdx.x] = 0;
			if (i < length)
				shared[threadIdx.x] += x[i] * y[i];
			if (i + blockDim.x < length)
				shared[threadIdx.x] += x[i + blockDim.x] * y[i + blockDim.x];

			__syncthreads();

			// Reduce pairs in shared memory.
			int start = pow2roundup((blockDim.x + blockDim.x % 2) >> 1);
			for (int s = start; s > warpSize; s >>= 1)
			{
				if ((threadIdx.x < s) && (threadIdx.x + s < half))
					shared[threadIdx.x] += shared[threadIdx.x + s];

				__syncthreads();
			}

			// Unroll last 32 iterations of loop (64 elements).
			// There is no need for synchronizations, since all accesses
			// are within single warp.
			if (threadIdx.x < warpSize)
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

		config_t configs[MAX_CONFIGS];
		double* buffer = NULL;
		size_t szbuffer = 0;
	}
}

// computes the inner product of x and y
// x and y are vectors of length N
inline double ss_dot(
	const double* const __restrict__ x, const double* const __restrict__ y, const int length)
{
	using namespace gpu;
	using namespace gpu::ss_dot_kernel;

	{
		size_t size = configs[0].grid.x + configs[0].grid.x % 2;
		if (szbuffer != size)
		{
			if (szbuffer)
				CUDA_ERR_CHECK(cudaFree(buffer));
			CUDA_ERR_CHECK(cudaMalloc(&buffer, sizeof(double) * size));
			szbuffer = size;
		}

		CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<
			configs[0].grid, configs[0].block, configs[0].block.x * sizeof(double)>>>(
			length, x, y, buffer));
	}

	for (int i = 1, szbuffer = configs[0].grid.x; szbuffer != 1; i++)
	{
		int x_length = szbuffer / 2 + szbuffer % 2;
		int y_length = szbuffer / 2;

		const double* x_dev = buffer;
		const double* y_dev = buffer + x_length;

		CUDA_LAUNCH_ERR_CHECK(ss_sum_kernel::kernel<1, gpu::double1><<<
			configs[i].grid, configs[i].block, configs[i].block.x * sizeof(double)>>>(
			x_length, x_dev, y_length, y_dev, buffer));
		
		szbuffer = configs[i].grid.x;
	}

	// record the number of floating point oporations
	flops_blas1 += 2 * length;

	CUDA_ERR_CHECK(cudaDeviceSynchronize());

	double result;
	CUDA_ERR_CHECK(cudaMemcpy(&result, &buffer[0], sizeof(double), cudaMemcpyDeviceToHost));
	return result;
}

namespace gpu
{
	namespace ss_norm2_kernel
	{
		// computes the 2-norm of x
		// x is a vector of length N
		template<short V, typename T>
		__global__ void kernel(const int length,
			const double* const __restrict__ x, double* const __restrict__ result)
		{
			extern __shared__ double shared[];

			int half = length / 2 + length % 2;

			// Each block hanldes (2 * blockDim.x) elements of reduction.
			int i = 2 * blockIdx.x * blockDim.x + threadIdx.x;

			// Load product of first 2 pairs into shared memory:
			// idx-th and (idx + blockDim.x)-th.
			shared[threadIdx.x] = 0;
			if (i < length)
				shared[threadIdx.x] += x[i] * x[i];
			if (i + blockDim.x < length)
				shared[threadIdx.x] += x[i + blockDim.x] * x[i + blockDim.x];

			__syncthreads();

			// Reduce pairs in shared memory.
			int start = pow2roundup((blockDim.x + blockDim.x % 2) >> 1);
			for (int s = start; s > warpSize; s >>= 1)
			{
				if ((threadIdx.x < s) && (threadIdx.x + s < half))
					shared[threadIdx.x] += shared[threadIdx.x + s];

				__syncthreads();
			}

			// Unroll last 32 iterations of loop (64 elements).
			// There is no need for synchronizations, since all accesses
			// are within single warp.
			if (threadIdx.x < warpSize)
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

		config_t configs[MAX_CONFIGS];
		double* buffer = NULL;
		size_t szbuffer = 0;
	}
}

// computes the 2-norm of x
// x is a vector of length N
inline double ss_norm2(const double* const __restrict__ x, const int length)
{
	using namespace gpu;
	using namespace gpu::ss_norm2_kernel;

	{
		size_t size = configs[0].grid.x + configs[0].grid.x % 2;
		if (szbuffer != size)
		{
			if (szbuffer)
				CUDA_ERR_CHECK(cudaFree(buffer));
			CUDA_ERR_CHECK(cudaMalloc(&buffer, sizeof(double) * size));
			szbuffer = size;
		}

		CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<
			configs[0].grid, configs[0].block, configs[0].block.x * sizeof(double)>>>(
			length, x, buffer));
	}

	for (int i = 1, szbuffer = configs[0].grid.x; szbuffer != 1; i++)
	{
		int x_length = szbuffer / 2 + szbuffer % 2;
		int y_length = szbuffer / 2;

		const double* x_dev = buffer;
		const double* y_dev = buffer + x_length;

		CUDA_LAUNCH_ERR_CHECK(ss_sum_kernel::kernel<1, gpu::double1><<<
			configs[i].grid, configs[i].block, configs[i].block.x * sizeof(double)>>>(
			x_length, x_dev, y_length, y_dev, buffer));
		
		szbuffer = configs[i].grid.x;
	}
	
	// record the number of floating point oporations
	flops_blas1 += 2 * length;

	CUDA_ERR_CHECK(cudaDeviceSynchronize());

	double result;
	CUDA_ERR_CHECK(cudaMemcpy(&result, &buffer[0], sizeof(double), cudaMemcpyDeviceToHost));
	return sqrt(result);
}

namespace gpu
{
	namespace ss_fill_kernel
	{
		// sets entries in a vector to value
		// x is a vector of length N
		// value is th
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ x, const double value, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T vv;
			for (int v = 0; v < V; v++)			
				vv.v[v] = value;
			((T*)x)[i] = vv;
		}

		config_t config;
	}
}

// sets entries in a vector to value
// x is a vector of length N
// value is th
inline void ss_fill(double* __restrict__ x, const double value, const int N)
{
	using namespace gpu;
	using namespace gpu::ss_fill_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(x, value, N));
}

////////////////////////////////////////////////////////////////////////////////
//  blas level 1 vector-vector operations
////////////////////////////////////////////////////////////////////////////////

namespace gpu
{
	namespace ss_axpy_kernel
	{
		// computes y := alpha*x + y
		// x and y are vectors of length N
		// alpha is a scalar
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ y, const double alpha,
			const double* const __restrict__ x, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T yy;
			T xx = T::ld(x, i);
			T zz = T::ld(y, i);
			for (int v = 0; v < V; v++)			
				yy.v[v] = alpha * xx.v[v] + zz.v[v];
			T::stcs(y, i, yy);
		}

		config_t config;
	}
}

// computes y := alpha*x + y
// x and y are vectors of length N
// alpha is a scalar
inline void ss_axpy(
	double* __restrict__ y, const double alpha, const double* const __restrict__ x, const int N)
{
	using namespace gpu;
	using namespace gpu::ss_axpy_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(y, alpha, x, N));

	// record the number of floating point oporations
	flops_blas1 += 2 * N;
}

namespace gpu
{
	namespace ss_add_scaled_diff_kernel
	{
		// computes y = x + alpha*(l-r)
		// y, x, l and r are vectors of length N
		// alpha is a scalar
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ y,
			const double* const __restrict__ x, const double alpha,
			const double* const __restrict__ l, const double* const __restrict__ r, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T yy;
			T xx = T::ld(x, i);
			T ll = T::ld(l, i);
			T rr = T::ld(r, i);
			for (int v = 0; v < V; v++)			
				yy.v[v] = xx.v[v] + alpha * (ll.v[v] - rr.v[v]);
			T::stcs(y, i, yy);
		}

		config_t config;
	}
}

// computes y = x + alpha*(l-r)
// y, x, l and r are vectors of length N
// alpha is a scalar
inline void ss_add_scaled_diff(double* __restrict__ y,
	const double* const __restrict__ x, const double alpha,
	const double* const __restrict__ l, const double* const __restrict__ r, const int N)
{
	using namespace gpu;
	using namespace ss_add_scaled_diff_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(y, x, alpha, l, r, N));

	// record the number of floating point oporations
	flops_blas1 += 3 * N;
}

namespace gpu
{
	namespace ss_scaled_diff_kernel
	{
		// computes y = alpha*(l-r)
		// y, l and r are vectors of length N
		// alpha is a scalar
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ y,
			const double alpha,	const double* const __restrict__ l,
			const double* const __restrict__ r, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T yy;
			T ll = T::ld(l, i);
			T rr = T::ld(r, i);
			for (int v = 0; v < V; v++)			
				yy.v[v] = alpha * (ll.v[v] - rr.v[v]);
			T::stcs(y, i, yy);
		}

		config_t config;
	}
}

// computes y = alpha*(l-r)
// y, l and r are vectors of length N
// alpha is a scalar
inline void ss_scaled_diff(double* __restrict__ y, const double alpha,
	const double* const __restrict__ l, const double* const __restrict__ r, const int N)
{
	using namespace gpu;
	using namespace gpu::ss_scaled_diff_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(y, alpha, l, r, N));

	// record the number of floating point oporations
	flops_blas1 += 2 * N;
}

namespace gpu
{
	namespace ss_scale_kernel
	{
		// computes y := alpha*x
		// alpha is scalar
		// y and x are vectors of length N
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ y, const double alpha,
			const double* const __restrict__ x, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T yy;
			T xx = T::ld(x, i);
			for (int v = 0; v < V; v++)			
				yy.v[v] = alpha * xx.v[v];
			T::stcs(y, i, yy);
		}

		config_t config;
	}
}

// computes y := alpha*x
// alpha is scalar
// y and x are vectors of length N
inline void ss_scale(double* __restrict__ y, const double alpha,
	const double* const __restrict__ x, const int N)
{
	using namespace gpu;
	using namespace gpu::ss_scale_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(y, alpha, x, N));

	// record the number of floating point oporations
	flops_blas1 += N;
}

namespace gpu
{
	namespace ss_lcomb_kernel
	{
		// computes linear combination of two vectors y := alpha*x + beta*z
		// alpha and beta are scalar
		// y, x and z are vectors of length N
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ y, const double alpha,
			const double* const __restrict__ x, const double beta,
			const double* const __restrict__ z, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T yy;
			T xx = T::ld(x, i);
			T zz = T::ld(z, i);
			for (int v = 0; v < V; v++)			
				yy.v[v] = alpha * xx.v[v] + beta * zz.v[v];
			T::stcs(y, i, yy);
		}

		config_t config;
	}
}

// computes linear combination of two vectors y := alpha*x + beta*z
// alpha and beta are scalar
// y, x and z are vectors of length N
inline void ss_lcomb(double* __restrict__ y, const double alpha,
	const double* const __restrict__ x, const double beta,
	const double* const __restrict__ z, const int N)
{
	using namespace gpu;
	using namespace gpu::ss_lcomb_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(y, alpha, x, beta, z, N));

	// record the number of floating point oporations
	flops_blas1 += 3 * N;
}

namespace gpu
{
	namespace ss_copy_kernel
	{
		// copy one vector into another y := x
		// x and y are vectors of length N
		template<short V, typename T>
		__global__ void kernel(double* __restrict__ y, const double* const __restrict__ x, const int N)
		{
			int i = blockDim.x * blockIdx.x + threadIdx.x;
			if (V * i >= N) return;

			T xx = T::ld(x, i);
			((T*)y)[i] = xx;
		}

		config_t config;
	}
}

// copy one vector into another y := x
// x and y are vectors of length N
inline void ss_copy(double* y, const double* const __restrict__ x, const int N)
{
	using namespace gpu;
	using namespace gpu::ss_copy_kernel;

	CUDA_LAUNCH_ERR_CHECK(kernel<2, gpu::double2><<<config.grid, config.block>>>(y, x, N));
}

// conjugate gradient solver
// solve the linear system A*x = b for x
// the matrix A is implicit in the objective function for the diffusion equation
// the value in x constitute the "first guess" at the solution
// x(N)
// ON ENTRY contains the initial guess for the solution
// ON EXIT  contains the solution
inline bool ss_cg(int N, double* __restrict__ x, const double* const __restrict__ b,
	const int maxiters, const double tol)
{
	if (!cg_initialized)
	{
		printf("INITIALIZING CG STATE\n");
		cg_init(N);
	}

	// epslion value use for matrix-vector approximation
	double eps	 = 1.e-8;
	double eps_inv = 1. / eps;

	// allocate memory for temporary storage
	ss_fill(Fx,	0.0, N);
	ss_fill(Fxold, 0.0, N);
	ss_copy(xold, x, N);

	// matrix vector multiplication is approximated with
	// A*v = 1/epsilon * ( F(x+epsilon*v) - F(x) )
	//	 = 1/epsilon * ( F(x+epsilon*v) - Fxold )
	// we compute Fxold at startup
	// we have to keep x so that we can compute the F(x+exps*v)
	diffusion(x, Fxold);

	// v = x + epsilon*x
	ss_scale(v, 1.0 + eps, x, N);

	// Fx = F(v)
	diffusion(v, Fx);

	// r = b - A*x
	// where A*x = (Fx-Fxold)/eps
	ss_add_scaled_diff(r, b, -eps_inv, Fx, Fxold, N);

	// p = r
	ss_copy(p, r, N);

	// rold = <r,r>
	double rold = ss_dot(r, r, N), rnew = rold;

	// check for convergence
	if (sqrt(rold) < tol)
		return true;

	int iter = 1; 
	for ( ; iter <= maxiters; iter++)
	{
		// Ap = A*p
		ss_lcomb(v, 1.0, xold, eps, p, N);
		diffusion(v, Fx);
		ss_scaled_diff(Ap, eps_inv, Fx, Fxold, N);

		// alpha = rold / p'*Ap
		double alpha = rold / ss_dot(p, Ap, N);

		// x += alpha*p
		ss_axpy(x, alpha, p, N);

		// r -= alpha*Ap
		ss_axpy(r, -alpha, Ap, N);

		// find new norm
		rnew = ss_dot(r, r, N);

		// test for convergence
		if (sqrt(rnew) < tol)
		{
			iters_cg += iter;
			return true;
		}

		// p = r + rnew.rold * p
		ss_lcomb(p, 1.0, r, rnew / rold, p, N);

		rold = rnew;
	}
	iters_cg += iter;

	printf("ERROR: CG failed to converge after %d iterations\n", maxiters);
	printf("	   achived tol = %E, required tol = %E\n", sqrt(rnew), tol);
	
	return false;
}

#endif // LINALG_H

