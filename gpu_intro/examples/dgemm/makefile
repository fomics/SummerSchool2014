default: dense.gpu.single dense.gpu.double

dense.gpu.single: dense.cpp
	g++ -DHAVE_SINGLE -g -O3 -I$(shell dirname $(shell which nvcc))/../include $< -o $@ -L$(shell dirname $(shell which nvcc))/../lib64 -L$(shell dirname $(shell which nvcc))/../lib -lcublas -lcudart -lrt

dense.gpu.double: dense.cpp
	g++ -DHAVE_DOUBLE -g -O3 -I$(shell dirname $(shell which nvcc))/../include $< -o $@ -L$(shell dirname $(shell which nvcc))/../lib64 -L$(shell dirname $(shell which nvcc))/../lib -lcublas -lcudart -lrt

dense.cpu.single: dense.cpp
	icpc -DHAVE_SINGLE -O3 -xHost -DCPU $< -o $@ -mkl -lrt

dense.cpu.double: dense.cpp
	icpc -DHAVE_DOUBLE -O3 -xHost -DCPU $< -o $@ -mkl -lrt

dense.mic.native.single: dense.cpp
	icpc -DHAVE_SINGLE -O3 -mmic -DCPU $< -o $@ -mkl -lrt

dense.mic.native.double: dense.cpp
	icpc -DHAVE_DOUBLE -O3 -mmic -DCPU $< -o $@ -mkl -lrt

clean:
	rm -rf dense.gpu.single dense.gpu.double dense.cpu.single dense.cpu.double dense.mic.native.single dense.mic.native.double

