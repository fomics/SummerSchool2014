all: compute_grid_3d_1 compute_grid_3d_2

compute_grid_3d_1: compute_grid_3d_1.cu
	nvcc -g -O3 -arch=sm_20 $< -o $@

compute_grid_3d_2: compute_grid_3d_2.cu
	nvcc -g -O3 -arch=sm_20 $< -o $@

clean:
	rm -rf compute_grid_3d_1 compute_grid_3d_2

