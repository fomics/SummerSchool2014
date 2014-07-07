all: cuda_gpu_test

cuda_gpu_test: cuda_gpu_test.cu
	nvcc -arch=sm_20 -O3 $< -o $@

clean:
	rm -rf cuda_gpu_test

