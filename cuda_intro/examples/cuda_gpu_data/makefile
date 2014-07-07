all: cuda_gpu_data

cuda_gpu_data: cuda_gpu_data.cu
	nvcc -O3 $< -o $@

clean:
	rm -rf cuda_gpu_data

