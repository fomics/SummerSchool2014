program dot

use omp_lib

implicit none

integer      :: N = 100000000
integer      :: i, num_threads
real(kind=8), dimension(:), allocatable  :: a, b
real(kind=8) :: time, sum, expected

allocate(a(N), b(N))

num_threads = omp_get_max_threads()
print *, 'dot of vectors with length ', N, ' with ', num_threads, ' threads'

do i=1,N
    a(i) = 1.0d0/2.0d0
    b(i) = i
enddo

time = -omp_get_wtime();

sum = 0.0d0
do i=1,n
    sum = sum + a(i) * b(i)
end do
time = time + omp_get_wtime()

expected = (N+1.0d0)*N/4.0d0;
print *, 'relative error ', abs(expected-sum)/expected
print *, 'took ', time, ' seconds'

deallocate(a, b)

end ! program dot

