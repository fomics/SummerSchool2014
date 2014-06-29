program pi_app

use omp_lib

implicit none

integer :: num_steps = 500000000
integer :: i
real(kind=8) :: x, pi, pi_reference, sum, w, time, error

pi=0.0
sum=0.0

print *, "using ", omp_get_max_threads(), " OpenMP threads"

! start timer
time = -omp_get_wtime()

w = 1.0d0/num_steps

do i=1,num_steps
    x = (i-0.5d0)*w;
    sum = sum + 4.0d0/(1.0d0+x*x);
enddo
time = time + omp_get_wtime()

! estimate pi
pi = w*sum

! calculate the relative error of our estimation
pi_reference = 3.141592653589793238462643383279502884d0
error = abs(pi-pi_reference)/pi_reference

print *, num_steps, " steps approximates pi as : ", pi, ", with relative error ", error
print *, "the solution took ", time, " seconds"

end ! program pi_app

