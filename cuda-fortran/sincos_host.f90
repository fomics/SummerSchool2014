subroutine cuda_err_check(istat)

use iso_c_binding

! Module where we define CUDA functions to use them from Fortran
use cudafor

implicit none

  integer, intent(in) :: istat

  ! Variables for error handling
  type(c_ptr) :: msg_ptr
  integer, parameter :: msg_len = 1024
  character(len=msg_len), pointer :: msg

  if (istat .ne. 0) then
      ! Get void* pointer to the error string
      msg_ptr = cudaGetErrorString(istat)
      ! Cast void* error string pointer to string
      call c_f_pointer(msg_ptr, msg)
      write(*, *) 'CUDA Error: ', &
          ! Take entire message string until first occurence of '\0' symbol
          ! -- end-of-string marker in C
          msg(1:index(msg, c_null_char))
      return
  endif

end subroutine cuda_err_check



subroutine sincos (nx, ny, nz)

use iso_c_binding

! Module where we define CUDA functions to use them from Fortran
use cudafor

implicit none

  ! Declare sincos kernel arguments with "target" -- for c_loc
  integer, target, intent(in) :: nx, ny, nz
  real, target, allocatable, dimension(:,:,:) :: h_x, h_y, h_xy1, h_xy2
  type(c_ptr) :: d_x, d_y, d_xy
    
  integer :: sum_host
  integer(c_int) :: istat

  integer, parameter :: szblock = 128

  ! CUDA constants
  integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
  integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2

  integer(c_size_t) :: szreal, szint, szptr, offset = 0

  type(dim3) :: blocks, threads

  integer :: i, j, k, step
  real, volatile :: start, finish

  szreal = sizeof(h_x(1, 1, 1))
  szint = sizeof(nx)
  szptr = sizeof(d_x)
  
  allocate(h_x(nx, ny, nz))
  allocate(h_y(nx, ny, nz))
  allocate(h_xy1(nx, ny, nz))
  allocate(h_xy2(nx, ny, nz))

  ! Fill input arrays with random values
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        call random_number(h_x(i,j,k))
        call random_number(h_y(i,j,k))
      enddo
    enddo
  enddo

  ! We repeat GPU logic twice to show that the first kernel invocation includes
  ! the time of GPU device initialization, which is done one time per
  ! application run.
  do step = 1, 2

    call cpu_time(start)

    ! Allocate memory on GPU.
    istat = cudaMalloc(d_x, nx * ny * nz * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(d_y, nx * ny * nz * szreal)
    call cuda_err_check(istat)
    istat = cudaMalloc(d_xy, nx * ny * nz * szreal)
    call cuda_err_check(istat)

	! XXX Several blocks of code below are equivalent to sum_kernel<<<...>>> in C.

    ! Create CUDA compute grid configuration.
    blocks = dim3(nx / szblock, ny, nz)
    if (mod(nx, szblock) .ne. 0) blocks%x = blocks%x + 1
    threads = dim3(szblock, 1, 1)

    ! Copy input data from CPU to GPU memory.
    istat = cudaMemcpy(d_x, c_loc(h_x), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)
    istat = cudaMemcpy(d_y, c_loc(h_y), nx * ny * nz * szreal, cudaMemcpyHostToDevice)
    call cuda_err_check(istat)

    ! Setup CUDA compute grid configuration.
    istat = cudaConfigureCall(blocks, threads, int8(0), c_null_ptr)
    call cuda_err_check(istat)

    ! Setup CUDA kernel arguments (individually)
    offset = 0
    istat = cudaSetupScalarArgument(c_loc(nx), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nx)
    istat = cudaSetupScalarArgument(c_loc(ny), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(ny)
    istat = cudaSetupScalarArgument(c_loc(nz), szint, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(nz)

    ! Previous arguments were 4 bytes in size, total offset is 12 bytes.
    ! Since the next 3 arguments are 8 bytes each and have to be aligned
    ! by 8 bytes boundary, it's necessary to insert a 4-byte spacing here
    ! (from 12 to 16).
    offset = offset + 4
    istat = cudaSetupArrayArgument(d_x, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(d_x)
    istat = cudaSetupArrayArgument(d_y, szptr, offset)
    call cuda_err_check(istat)
    offset = offset + sizeof(d_y)
    istat = cudaSetupArrayArgument(d_xy, szptr, offset)
    call cuda_err_check(istat)

    ! Finally, launch CUDA kernel after we configured everything.
    ! Kernel is identified by C-string with its name (must be
    ! undecorated, i.e. with extern "C")
    istat = cudaLaunch('sincos_kernel' // c_null_char)
    call cuda_err_check(istat)

    ! Copy back results from GPU to CPU memory.
    istat = cudaMemcpy(c_loc(h_xy1), d_xy, nx * ny * nz * szreal, cudaMemcpyDeviceToHost)
    call cuda_err_check(istat)

    ! Free allocated GPU memory.
    istat = cudaFree(d_x)
    call cuda_err_check(istat)
    istat = cudaFree(d_y)
    call cuda_err_check(istat)
    istat = cudaFree(d_xy)
    call cuda_err_check(istat)

    call cpu_time(finish)

    print *, 'gpu time = ', finish - start

  enddo

  call cpu_time(start)

  ! Control CPU implementation
  !$omp parallel do
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        h_xy2(i,j,k) = sin(h_x(i,j,k)) + cos(h_y(i,j,k))
      enddo
    enddo
  enddo

  call cpu_time(finish)

  print *, 'cpu time = ', finish - start

  ! Compare results
  print *, 'diff = ', maxval(h_xy1 - h_xy2)

  deallocate(h_x)
  deallocate(h_y)
  deallocate(h_xy1)
  deallocate(h_xy2)

end subroutine sincos

