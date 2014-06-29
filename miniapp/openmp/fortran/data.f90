module data

implicit none

! define some helper types that can be used to pass simulation
! data around without haveing to pass individual parameters

type discretizationT
    integer nx      ! x dimension
    integer ny      ! y dimension
    integer nt      ! number of time steps
    integer N       ! total number of grid points
    real (kind=8) dt    ! time step size
    real (kind=8) dx    ! distance between grid points
    real (kind=8) alpha ! dx^2/(D*dt)
end type

! store information about the domain decomposition
type subdomainT
    ! i and j dimensions of the global decomposition
    integer i_dimension
    integer j_dimension
    ! the i and j index of this sub-domain
    integer i_index
    integer j_index
    ! boolean flags that indicate whether the sub-domain is on 
    ! any of the 4 global boundaries
    logical on_boundary_north
    logical on_boundary_south
    logical on_boundary_east
    logical on_boundary_west
end type

! fields that hold the solution
real (kind=8), allocatable :: x_new(:,:), x_old(:,:)
real (kind=8), allocatable :: bndN(:), bndE(:), bndS(:), bndW(:)
type(discretizationT) :: options

end module data
