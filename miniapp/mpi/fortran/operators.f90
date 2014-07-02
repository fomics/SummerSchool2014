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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do the boundary exchange
    ! TODOSS
    ! * uncomment to get the first 1/2 of communications working (don't forget about the MPI_Waitall below!)
    ! * try visualizing: What is missing?
    ! * finish communications
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    num_requests = 0
    ! if (domain%neighbour_north>=0) then
    !     ! set tag to be the sender's rank
    !     ! post receive
    !     call mpi_irecv(bndN, nx, MPI_DOUBLE, domain%neighbour_north, domain%neighbour_north, &
    !         MPI_COMM_WORLD, requests(num_requests+1), err)

    !     ! pack north buffer
    !     buffN = u(:,ny)

    !     ! post send
    !     call mpi_isend(buffN, nx, MPI_DOUBLE, domain%neighbour_north, domain%rank, &
    !         MPI_COMM_WORLD, requests(num_requests+2), err)

    !     num_requests = num_requests + 2
    ! endif
    ! if (domain%neighbour_south>=0) then
    !     ! post receive
    !     call mpi_irecv(bndS, nx, MPI_DOUBLE, domain%neighbour_south, domain%neighbour_south, &
    !         MPI_COMM_WORLD, requests(num_requests+1), err)

    !     ! pack north buffer
    !     buffS = u(:,1)

    !     ! post send
    !     call mpi_isend(buffS, nx, MPI_DOUBLE, domain%neighbour_south, domain%rank, &
    !         MPI_COMM_WORLD, requests(num_requests+2), err)

    !     num_requests = num_requests + 2
    ! endif
    
    ! wait on the receives
    ! TODOSS
    ! * fill in the blanks
    ! * can this wait be moved?
    ! call mpi_waitall( .. , .. , .. , .. )
    
    ! the interior grid points
    do j = 2, jend
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
    do j = 2, jend
        s(i,j) = -(4.+alpha) * u(i,j)        &
                    + u(i-1, j) + u(i, j-1) + u(i, j+1) &
                    + alpha*x_old(i,j) + bndE(j) &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    ! the west boundary
    i = 1
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
    do i = 2, iend
        s(i,j) = -(4.+alpha) * u(i,j)        &
                    + u(i-1, j) + u(i+1, j) + u(i, j-1) &
                    + alpha*x_old(i,j) + bndN(i) &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    i = options%nx ! NE corner
    s(i,j) = -(4.+alpha) * u(i,j)       &
                + u(i-1, j) + u(i, j-1) &
                + alpha*x_old(i,j)      &
                + bndE(j) + bndN(i) &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    ! the south boundary
    j = 1
    i = 1 ! SW corner
    s(i,j) = -(4.+alpha) * u(i,j)       &
                + u(i+1, j) + u(i, j+1) &
                + alpha*x_old(i,j)      &
                + bndW(j) + bndS(i) &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    ! south boundary
    do i = 2, iend
        s(i,j) = -(4.+alpha) * u(i,j)           &
                    + u(i-1,j  ) + u(i+1,j  )   &
                                 + u(i  ,j+1)   &
                    + alpha*x_old(i,j)          &
                    + bndS(i) &
                    + dxs*u(i,j)*(1.0_8 - u(i,j))
    end do

    i = options%nx ! SE corner
    s(i,j) = -(4.+alpha) * u(i,j)       &
                + u(i-1,j) + u(i, j+1)  &
                + alpha*x_old(i,j)      &
                + bndE(j) + bndS(i) &
                + dxs*u(i,j)*(1.0_8 - u(i,j))

    ! accumulate the flop counts
    ! 8 ops total per point
    flops_diff =  flops_diff                    &
                    + 12 * (options%nx-2) * (options%ny-2) & ! interior points
                    + 11 * (options%nx-2 + options%ny-2) &   ! NESW boundary points
                    + 11 * 4                                 ! corner points
end

end module operators

