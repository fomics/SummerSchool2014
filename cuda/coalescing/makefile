all: coalescing coalescing_res coalescing_comp

coalescing.o: coalescing.cu
	nvcc -arch=sm_20 -c $< -o $@

coalescing: coalescing.o 
	nvcc -arch=sm_20 $< -o $@

coalescing_res.o: coalescing_res.cu
	nvcc -arch=sm_20 -c $< -o $@

coalescing_res: coalescing_res.o 
	nvcc -arch=sm_20 $< -o $@

coalescing_comp.o: coalescing_comp.cu
	nvcc -arch=sm_20 -c $< -o $@

coalescing_comp: coalescing_comp.o 
	nvcc -arch=sm_20 $< -o $@

clean:
	rm -rf *.o coalescing coalescing_res coalescing_comp
