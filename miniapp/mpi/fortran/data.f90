module data

implicit none

! define some helper types that can be used to pass simulation
! data around without haveing to pass individual parameters

type discretizationT
    integer nx          ! x dimension
    integer ny          ! y dimension
    integer global_nx   ! x dimension of global domain
    integer global_ny   ! y dimension of global domain
    integer nt          ! number of time steps
    integer N           ! total number of grid points
    real (kind=8) dt    ! time step size
    real (kind=8) dx    ! distance between grid points
    real (kind=8) alpha ! dx^2/(D*dt)
end type

! store information about the domain decomposition
type subdomainT
    ! i and j dimensions of the global decomposition
    integer ndomx
    integer ndomy

    ! the i and j index of this sub-domain
    integer domx
    integer domy

    ! the i and j bounding box of this sub-domain
    integer startx
    integer starty
    integer endx
    integer endy

    ! the rank of neighbouring domains
    integer neighbour_north
    integer neighbour_east
    integer neighbour_south
    integer neighbour_west

    ! mpi info
    integer size
    integer rank
end type

! fields that hold the solution
real (kind=8), allocatable :: x_new(:,:), x_old(:,:)
real (kind=8), allocatable :: bndN(:), bndE(:), bndS(:), bndW(:)
real (kind=8), allocatable :: buffN(:), buffS(:), buffE(:), buffW(:)

type(discretizationT) :: options
type(subdomainT)      :: domain

end module data
