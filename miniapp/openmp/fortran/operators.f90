!******************************************
! operators.f90
! based on min-app code written by Oliver Fuhrer, MeteoSwiss
! modified by Ben Cumming, CSCS
! *****************************************

! Description: Contains simple operators which can be used on 3d-meshes

module operators

use stats,     only: flops_diff
use data,      only: discretizationT, x_old, options, bndN, bndE, bndS, bndW

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
    integer :: iend, jend

    dxs   = 1000.*(options%dx ** 2)
    alpha = options%alpha
    iend  = options%nx-1
    jend  = options%ny-1

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
end subroutine

end module operators

