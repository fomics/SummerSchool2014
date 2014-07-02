#!/bin/bash

OMP_NUM_THREADS=8 aprun -n 1 -cc none ./main 256 256 50 0.001 | grep second
OMP_NUM_THREADS=4 aprun -n 2 -cc none ./main 256 256 50 0.001 | grep second
OMP_NUM_THREADS=2 aprun -n 4 -cc none ./main 256 256 50 0.001 | grep second
OMP_NUM_THREADS=1 aprun -n 8          ./main 256 256 50 0.001 | grep second
