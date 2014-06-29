program sum_threads

use omp_lib

implicit none

    integer :: num_threads, expected
    integer :: sum=0

    num_threads = omp_get_max_threads()
    print *, 'sum with ', num_threads, ' threads'

    !$omp parallel
    sum = sum + omp_get_thread_num()+1
    !$omp end parallel

    ! use formula for sum of arithmetic sequence: sum(1:n) = (n+1)*n/2
    expected = (num_threads+1)*num_threads/2
    if (sum==expected) then
        print *, "sum ", sum, ' matches the expected value'
    else
        print *, "sum ", sum, ' does not match the expected value ', expected
    endif

end ! program sum_threads

