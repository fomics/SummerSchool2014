all: dot.exe normalize.exe hello_world.exe sum_threads.exe pi.exe

pi.exe : pi.f90
	ftn pi.f90 -fopenmp -o pi.exe

normalize.exe : normalize.f90
	ftn normalize.f90 -fopenmp -o normalize.exe

hello_world.exe : hello_world.f90
	ftn hello_world.f90 -fopenmp -o hello_world.exe

dot.exe : dot.f90
	ftn dot.f90 -fopenmp -o dot.exe

sum_threads.exe : sum_threads.f90
	ftn sum_threads.f90 -fopenmp -o sum_threads.exe

clean :
	rm -f *.exe

