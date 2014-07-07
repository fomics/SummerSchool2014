!******************************************
! operators.f90
! based on min-app code written by Oliver Fuhrer, MeteoSwiss
! modified by Ben Cumming, CSCS
! *****************************************

! Description: Contains simple operators which can be used on 3d-meshes

module operators

use mpi
use stats,     only: flops_diff
use data,      only: discretizationT, x_old, options, bndN, bndE, bndS, bndW, domain, buffN, buffS, buffE, buffW

implicit none

contains

!==============================================================================

subroutine diffusion(u, s)
    ! arguments
    real (kind=8), intent(in)  :: u(options%nx, options%ny)
    real (kind=8), intent(out) :: s(options%nx, options%ny)

    ! local variables
    real (kind=8) :: alpha, dxs
    integer :: i, j
    integer :: iend, jend, nx, ny
    integer :: stats(MPI_STATUS_SIZE,8)
    integer :: requests(8)
    integer :: num_requests, err

    dxs   = 1000.*(options%dx ** 2)
    alpha = options%alpha
    iend  = options%nx-1
    jend  = options%ny-1
    nx  = options%nx
    ny  = options%ny

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$acc data present(u, buffN, buffS, buffE, buffW)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !$acc parallel
    !$acc loop
    do i = 1,nx
        buffN(i) = u(i,ny)
    end do
    !$acc loop
    do i = 1,nx
        buffS(i) = u(i,1)
    end do
    !$acc loop
    do j = 1,ny
        buffE(j) = u(nx,j)
    end do
    !$acc loop
    do j = 1,ny
        buffW(j) = u(1,j)
    end do
    !$acc end parallel

#ifdef USE_G2G
    !$acc host_data use_device(bndN, buffN, bndS, buffS, bndE, buffE, bndW, buffW)
#else
    !$acc update host(buffN, buffE, buffS, buffW)
#endif
    num_requests = 0
    if (domain%neighbour_north>=0) then
        call mpi_irecv(bndN, nx, MPI_DOUBLE, domain%neighbour_north, domain%neighbour_north, &
            MPI_COMM_WORLD, requests(num_requests+1), err)
        call mpi_isend(buffN, nx, MPI_DOUBLE, domain%neighbour_north, domain%rank, &
            MPI_COMM_WORLD, requests(num_requests+2), err)

        num_requests = num_requests + 2
    endif
    if (domain%neighbour_south>=0) then
        call mpi_irecv(bndS, nx, MPI_DOUBLE, domain%neighbour_south, domain%neighbour_south, &
            MPI_COMM_WORLD, requests(num_requests+1), err)
        call mpi_isend(buffS, nx, MPI_DOUBLE, domain%neighbour_south, domain%rank, &
            MPI_COMM_WORLD, requests(num_requests+2), err)

        num_requests = num_requests + 2
    endif
    if (domain%neighbour_east>=0) then
        call mpi_irecv(bndE, ny, MPI_DOUBLE, domain%neighbour_east, domain%neighbour_east, &
            MPI_COMM_WORLD, requests(num_requests+1), err)
        call mpi_isend(buffE, ny, MPI_DOUBLE, domain%neighbour_east, domain%rank, &
            MPI_COMM_WORLD, requests(num_requests+2), err)

        num_requests = num_requests + 2
    endif
    if (domain%neighbour_west>=0) then
        call mpi_irecv(bndW, ny, MPI_DOUBLE, domain%neighbour_west, domain%neighbour_west, &
            MPI_COMM_WORLD, requests(num_requests+1), err)
        ! post send
        call mpi_isend(buffW, ny, MPI_DOUBLE, domain%neighbour_west, domain%rank, &
            MPI_COMM_WORLD, requests(num_requests+2), err)

        num_requests = num_requests + 2
    endif
#ifdef USE_G2G
    !$acc end host_data
#endif
    call mpi_waitall(num_requests, requests, stats, err)
#ifndef USE_G2G
    !$acc update device (bndE, bndW, bndS, bndN)
#endif

    !$acc parallel 
    ! the interior grid points
    !$acc loop
    do j = 2, jend
        !$acc loop
        do i = 2, iend
            s(i,j) = -(4.+alpha) * u(i,j)           &   ! central point
                        + u(i-1, j) + u(i+1, j)     &   ! east and west
                        + u(i, j-1) + u(i, j+1)     &   ! north and south
                        + alpha*x_old(i,j) &
                        + dxs*u(i,j)*(1.0_8 - u(i,j))
        end do
    end do

    ! the east boundary
    i = options%nx
    !$acc loop
    do j = 2, jend
        s(i,j) = -(4.+alpha) * u(i,j)        &
                    + u(i-1, j) + u(i, j-1) + u(i, j+1) &
                    + alpha*x_old(i,j) + bndE(j) &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    ! the west boundary
    i = 1
    !$acc loop
    do j = 2, jend
        s(i,j) = -(4.+alpha) * u(i,j)         &
                    + u(i+1, j) + u(i, j-1) + u(i, j+1) &
                    + alpha*x_old(i,j) + bndW(j) &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    ! the north boundary (plus NE and NW corners)
    j = options%ny
    i = 1 ! NW corner
    s(i,j) = -(4.+alpha) * u(i,j)           &
                + u(i+1, j) + u(i, j-1)     &
                + alpha*x_old(i,j)          &
                + bndW(j) + bndN(i) &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    ! north boundary
    !$acc loop
    do i = 2, iend
        s(i,j) = -(4.+alpha) * u(i,j)                   &
                    + u(i-1, j) + u(i+1, j) + u(i, j-1) &
                    + alpha*x_old(i,j) + bndN(i)        &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    i = options%nx ! NE corner
    s(i,j) = -(4.+alpha) * u(i,j)       &
                + u(i-1, j) + u(i, j-1) &
                + alpha*x_old(i,j)      &
                + bndE(j) + bndN(i)     &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    ! the south boundary
    j = 1
    i = 1 ! SW corner
    s(i,j) = -(4.+alpha) * u(i,j)       &
                + u(i+1, j) + u(i, j+1) &
                + alpha*x_old(i,j)      &
                + bndW(j) + bndS(i)     &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    ! south boundary
    !$acc loop
    do i = 2, iend
        s(i,j) = -(4.+alpha) * u(i,j)           &
                    + u(i-1,j  ) + u(i+1,j  )   &
                                 + u(i  ,j+1)   &
                    + alpha*x_old(i,j)          &
                    + bndS(i)                   &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    i = options%nx ! SE corner
    s(i,j) = -(4.+alpha) * u(i,j)       &
                + u(i-1,j) + u(i, j+1)  &
                + alpha*x_old(i,j)      &
                + bndE(j) + bndS(i)     &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$acc end parallel
    !$acc end data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! accumulate the flop counts
    ! 8 ops total per point
    flops_diff =  flops_diff                                &
                    + 12 * (options%nx-2) * (options%ny-2)  &! interior points
                    + 11 * (options%nx-2 + options%ny-2)    &! NESW boundary points
                    + 11 * 4                                 ! corner points
end

end module operators

