all: reduction

reduction: reduction.cu
	nvcc -g -O3 -arch=sm_35 $< -o $@

clean:
	rm -rf reduction

