// *****************************************
// operators.f90
// based on min-app code written by Oliver Fuhrer, MeteoSwiss
// modified by Ben Cumming, CSCS
// *****************************************

// Description: Contains simple operators which can be used on 3d-meshes

#ifndef OPERATORS_H
#define OPERATORS_H

#include "check.h"
#include "data.h"
#include "stats.h"

#define U(j,i)    up[(i) + (j)*nx]
#define S(j,i)    sp[(i) + (j)*nx]
#define X(j,i) x_old[(i) + (j)*nx]

// TODO
// Taget: Given that up and sp data are already on GPU, implement a set of
// GPU kernels for diffusion operator
// Things to note:
// 1) Look into how kernels in linalg are implemented
// 2) Note the code which provides kernels with optimal compute grid configs
// 3) Note options structure has cpu:: and gpu:: versions
// 4) Note you can turn every loop of diffusion into separate GPU kernel
// and run all of them asynchronously.
namespace gpu
{
	namespace diffusion_interior_grid_points_kernel
	{
		template<short V, typename T>
		__global__ void kernel(const double* const __restrict__ up, double* __restrict__ sp)
		{
			// TODO implement GPU kernel, equivalent to the original code in cuda-1/operators.c
		}

		config_t config;
	}

	namespace diffusion_east_west_boundary_points_kernel
	{
		template<short V, typename T>
		__global__ void kernel(const double* const __restrict__ up, double* __restrict__ sp)
		{
			// TODO implement GPU kernel, equivalent to the original code in cuda-1/operators.c
		}

		config_t config;
	}

	namespace diffusion_north_south_boundary_points_kernel
	{
		template<short V, typename T>
		__global__ void kernel(const double* const __restrict__ up, double* __restrict__ sp)
		{
			// TODO implement GPU kernel, equivalent to the original code in cuda-1/operators.c
		}

		config_t config;
	}

	namespace diffusion_corner_points_kernel
	{
		__global__ void kernel(const double* const __restrict__ up, double* __restrict__ sp)
		{
			// TODO implement GPU kernel, equivalent to the original code in cuda-1/operators.c
		}
	}
}

inline void diffusion(const double* const __restrict__ up, double* __restrict__ sp)
{
	{
		using namespace gpu;

		// Launch kernel for parallel processing of interior points.
		{
			using namespace diffusion_interior_grid_points_kernel;
			CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<config.grid, config.block>>>(up, sp));
		}

		// Launch kernels for parallel processing of boundary points.
		{
			using namespace diffusion_east_west_boundary_points_kernel;
			CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<config.grid, config.block>>>(up, sp));
		}
		{
			using namespace diffusion_north_south_boundary_points_kernel;
			CUDA_LAUNCH_ERR_CHECK(kernel<1, gpu::double1><<<config.grid, config.block>>>(up, sp));
		}

		// Finally, single-threaded processing of corner points.
		{
			using namespace diffusion_corner_points_kernel;
			CUDA_LAUNCH_ERR_CHECK(kernel<<<1, 1>>>(up, sp));
		}
	}

	{
		using namespace cpu;
		
		// Accumulate the flop counts
		// 8 ops total per point
		flops_diff +=
			+ 12 * (options.nx - 2) * (options.ny - 2) // interior points
			+ 11 * (options.nx - 2  +  options.ny - 2) // NESW boundary points
			+ 11 * 4;                                  // corner points}
	}
}

#endif // OPERATORS_H

