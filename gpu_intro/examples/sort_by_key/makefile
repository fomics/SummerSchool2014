#TBBINC =
#TBBLIB = /opt/intel/tbb/lib/intel64/gcc4.4
#INTELLIB = /opt/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64/

#TBBINC = /opt/intel/14.0.1.106/composer_xe_2013_sp1.1.106/tbb/include
#TBBLIB = /opt/intel/14.0.1.106/composer_xe_2013_sp1.1.106/tbb/lib/intel64/gcc4.4
#INTELLIB = /opt/intel/14.0.1.106/composer_xe_2013_sp1.1.106/compiler/lib/intel64/

TBBINC = /apps/dommic/intel/tbb/include
TBBLIB = /apps/dommic/intel/tbb/lib/intel64/
INTELLIB = /apps/dommic/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64

all: \
	sbk_soa_thrust.cpu sbk_aos_thrust.cpu \
	sbk_soa_thrust.gpu sbk_aos_thrust.gpu \
	sbk_aos_tbb.cpu sbk_aos_tbb.mic

sbk_soa_thrust.cpu: sbk_soa_thrust.cu
	nvcc -O3 -DCPU -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP $< -o $@ -lrt

sbk_aos_thrust.cpu: sbk_aos_thrust.cu
	nvcc -O3 -DCPU -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP $< -o $@ -lrt

sbk_soa_thrust.gpu: sbk_soa_thrust.cu
	nvcc -arch=sm_35 -O3 $< -o $@ -lrt

sbk_aos_thrust.gpu: sbk_aos_thrust.cu
	nvcc -arch=sm_35 -O3 $< -o $@ -lrt

sbk_aos_tbb.cpu: sbk_aos_tbb.cpp
	icpc -std=c++11 -I$(TBBINC) -ipo -O3 -no-prec-div -xHost -DCPU $< -o $@ -lrt -L$(TBBLIB) -ltbb -Wl,-rpath=$(TBBLIB) -Wl,-rpath=$(INTELLIB)

sbk_aos_tbb.mic: sbk_aos_tbb.cpp
	icpc -std=c++11 -wd3218 -wd2568 -I$(TBBINC) -ipo -O3 -no-prec-div -xHost $< -o $@ -lrt -L$(TBBLIB) -tbb -Wl,-rpath=$(TBBLIB) -Wl,-rpath=$(INTELLIB)

sbk_aos_tbb.mic.native: sbk_aos_tbb.cpp
	icpc -std=c++11 -I$(TBBINC) -mmic -openmp -O3 -DCPU $< -o $@ -lrt -L$(TBBLIB)/../mic -tbb -Wl,-rpath=$(TBBLIB)/../mic -Wl,-rpath=$(INTELLIB)/../mic

clean:
	rm -rf \
	sbk_soa_thrust.cpu sbk_aos_thrust.cpu \
	sbk_soa_thrust.gpu sbk_aos_thrust.gpu \
	sbk_aos_tbb.cpu sbk_aos_tbb.mic

test_mic:
	MIC_LD_LIBRARY_PATH=/users/mikushin/forge/kernelgen/doc/hpcac_swiss_2014/gpu_intro/sdk/sort_by_key/:/apps/dommic/intel/composer_xe_2013.2.146/compiler/lib/mic ./benchmark ./sbk_aos_tbb.mic 1024 $$((1024*1024*256)) 2 >results/sbk_aos_tbb.mic.xeon-phi-5110p
