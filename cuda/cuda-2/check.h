#ifndef CHECK_H
#define CHECK_H

#if 1
#define CUDA_ERR_CHECK(x) x;
#define CUDA_LAUNCH_ERR_CHECK(...) __VA_ARGS__;
#else
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
#endif

#endif // CHECK_H

