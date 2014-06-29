all: dot.exe normalize.exe hello_world.exe sum_threads.exe pi.exe

pi.exe : pi.cpp
	CC pi.cpp -fopenmp -o pi.exe

normalize.exe : normalize.cpp
	CC normalize.cpp -fopenmp -o normalize.exe

hello_world.exe : hello_world.cpp
	CC hello_world.cpp -fopenmp -o hello_world.exe

dot.exe : dot.cpp
	CC dot.cpp -fopenmp -o dot.exe

sum_threads.exe : sum_threads.cpp
	CC sum_threads.cpp -fopenmp -o sum_threads.exe

clean :
	rm -f *.exe

