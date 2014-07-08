#include "data.h"

#include <stdio.h>

struct discretization_t cpu::options;

cudaDeviceProp cpu::props;

// fields that hold the solution
__constant__ double *gpu::x_old = NULL; // 2d
__constant__ double *gpu::bndN = NULL, *gpu::bndE = NULL;
__constant__ double *gpu::bndS = NULL, *gpu::bndW = NULL;

__constant__ struct discretization_t gpu::options;

