!******************************************
! implicit time stepping implementation of 2D diffusion problem
! Ben Cumming, CSCS
! *****************************************

! A small benchmark app that solves the 2D fisher equation using second-order
! finite differences.

! Syntax:
!   ./main nx ny nt

program diffusion_serial

    ! modules
    use mpi
    use omp_lib
    use io,     only: write_header, write_parallel
    use stats,  only: flops_diff, flops_bc, flops_blas1, iters_cg, iters_newton
    use linalg, only: ss_copy, ss_scale, ss_cg, ss_axpy, ss_norm2, cg_init, r, Ap, p, Fx, Fxold, v, xold
    use data,   only: subdomainT, discretizationT, x_new, x_old, bndN, bndE, bndS, bndW, options, domain, buffN, buffW, buffE, buffS
    use operators,    only: diffusion

    implicit none

    ! variables
    integer :: nx, ny, nt, N

    real (kind=8), allocatable :: b(:), deltax(:)
    real (kind=8) :: timespent, time_in_bcs, time_in_diff, time_in_other
    real (kind=8) :: alpha, residual
    real (kind=8) :: one, zero
    real (kind=8) :: tolerance
    real (kind=8) :: x, y, xc, yc, radius

    integer :: ierr, timestep
    integer :: i, j, it
    integer :: istart, iend, jstart, jend
    integer :: iseed(80), nseed
    integer :: flops_total
    integer :: output
    integer :: err

    logical :: converged, cg_converged
    logical :: verbose_output

    ! ****************** read command line arguments ******************

    ! read command line arguments
    call readcmdline(options)

    call initialize_mpi(options, domain)

    nx = options%nx
    ny = options%ny
    N  = options%N
    nt = options%nt

    ! ****************** setup compute domain ******************

    if (domain%rank == 0) then
        write(*,'(A)') '========================================================================'
        print *,       '                      Welcome to mini-stencil!'
        print *, 'mesh :: ', nx, '*', ny, '    dx =', options%dx
        print *, 'time :: ', nt, 'time steps from 0 .. ', options%nt*options%dt
        write(*,'(A)') '========================================================================'
    endif

    ! ****************** constants ******************
    one = 1.
    zero = 0.

    ! ****************** allocate memory ******************

    ! allocate global fields
    allocate(x_new(nx,ny), x_old(nx,ny), b(N), deltax(N), stat=ierr)
    call error(ierr /= 0, 'Problem allocating memory')
    allocate(bndN(nx), bndS(nx), stat=ierr)
    call error(ierr /= 0, 'Problem allocating memory')
    allocate(bndE(ny), bndW(ny), stat=ierr)
    call error(ierr /= 0, 'Problem allocating memory')

    ! allocate buffers for MPI
    allocate(buffN(nx), buffS(nx), buffE(ny), buffW(ny))

    ! initialise memory required by conjugate gradient
    ! do this here so that is available outside the main OpenACC data region
    call cg_init(N)

    ! ****************** initialization ******************
    ! set dirichlet boundary conditions to 0 all around
    bndN  = 0
    bndS  = 0
    bndE  = 0
    bndW  = 0
    ! north west corner
    !bndN(1) = 0.2
    !bndW(ny) = 0.2
    ! south west corner
    !bndS(1) = 1
    !bndW(1) = 1
    ! north east corner
    !bndN(nx) = 0.1
    !bndE(ny) = 0.1
    ! south east corner
    !bndS(nx) = 0.5
    !bndE(1)  = 0.5

    ! set the initial condition
    ! a circle of concentration 0.1 centred at (xdim/4, ydim/4) with radius
    ! no larger than 1/8 of both xdim and ydim
    x_new = 0.0
    xc = 1.0/4.0
    yc = real(options%global_ny-1)*options%dx / 4
    radius = min(xc,yc)/2.0
    do j = domain%starty, domain%endy
        y = real(j-1)*options%dx
        do i = domain%startx, domain%endx
            x = real(i-1)*options%dx
            if ( (x-xc)**2 + (y-yc)**2 < radius**2) then
                x_new(i-domain%startx+1, j-domain%starty+1) = 0.1
            endif
        enddo
    enddo

    ! ****************** mpi reference version ******************
    time_in_bcs  = 0.0
    time_in_diff = 0.0
    flops_bc     = 0
    flops_diff   = 0
    flops_blas1  = 0
    iters_cg     = 0
    iters_newton = 0


    ! main timeloop
    alpha = options%alpha
    tolerance = 1.e-6

    converged = .true.

    timestep = 1;
    !$acc data copy(x_old, x_new, deltax, Ap, p, r, b, v, Fx, Fxold, x, xold, bndN, bndE, bndS, bndW, options, buffN, buffE, buffW, buffS)
    !$acc wait

    ! start timer
    timespent = -omp_get_wtime();

    do while(timestep<=nt .and. converged)
        ! set x_new and x_old to be the solution
        call ss_copy(x_old, x_new, N)

        converged = .false.
        cg_converged = .true.

        it = 0
        ! a slightly convoluted control statement on our while loop... but this
        ! is required to satisfy OpenACC's requirement that there should be not
        ! break/goto/return statements that leave a data region
        do while( (it<50) .and. (.not. converged) .and. cg_converged )
            ! compute residual : requires both x_new and x_old
            call diffusion(x_new, b)
            residual = ss_norm2(b, N)

            ! check for convergence
            if(residual<tolerance) then
                converged = .true.
            else    ! no convergence, so solve linear system
                ! solve linear system to get -deltax
                call ss_cg(deltax, b, 200, tolerance, cg_converged)

                ! update solution
                call ss_axpy(x_new, -one, deltax, N)
            endif
            it = it+1
        end do
        iters_newton = iters_newton + it

        ! output some statistics
        if (domain%rank==0) then
            if (converged .and. verbose_output) then
                write(*,*) 'step ', timestep, &
                        ' required ', it,  &
                        ' iterations for residual', residual
            endif
            if (.not. converged) then
                write(*,*) '!!!!!!!!!!!!!!!! step ', timestep, &
                        ' ERROR : nonlinear iterations failed to converge'
            endif
        endif

        timestep = timestep+1
    end do
    !$acc wait
    ! get times
    timespent = timespent + omp_get_wtime();
    !$acc end data

    flops_total = flops_diff + flops_blas1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write final solution to BOV file for visualization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! binary data
    call write_header('output.bov')
    call write_parallel('output.bin', x_new)

    ! print table sumarizing results
    if (domain%rank == 0) then
        output=20
        write(*,'(A)') '--------------------------------------------------------------------------------'
        write(*,*) 'simulation took ', timespent , ' seconds'
        write(*,*) iters_cg , ' conjugate gradient iterations', iters_cg/timespent, ' per second'
        write(*,*) iters_newton , ' nonlinear newton iterations'
        write(*,'(A)') '-------------------------------------------------------------------------------'
    end if


    ! ****************** cleanup ******************

    ! deallocate global fields
    deallocate(x_new, x_old)
    deallocate(bndN, bndS, bndE, bndW)

    call mpi_finalize(err)

contains

!==============================================================================

! swap data_in and data_out fields
! note: this is achieved without copying by using the Fortran move_alloc intrinsic
subroutine swap_data()
    implicit none

    real (kind=8), allocatable :: tmp(:,:)

    call move_alloc(x_new, tmp)
    call move_alloc(x_old, x_new)
    call move_alloc(tmp, x_old)
end

!==============================================================================

! read command line arguments
subroutine readcmdline(options)
    implicit none

    ! arguments
    type(discretizationT), intent(out) :: options

    ! local
    character(len=256) :: sarg
    integer :: nx, ny, nt
    integer :: nargs
    real (kind=8) :: t

    nargs = command_argument_count()
    if ( nargs < 4 .or. nargs > 5) then
        if(domain%rank == 0) then
            write(*,*) 'Usage: main nx ny nt t verbose'
            write(*,*) '  nx  number of gridpoints in x-direction'
            write(*,*) '  ny  number of gridpoints in y-direction'
            write(*,*) '  nt  number of timesteps'
            write(*,*) '  t   total time'
            write(*,*) '  verbose   (optional) if set verbose output is enabled'
        endif
        stop
    end if

    ! read nx
    call get_command_argument(1, sarg)
    nx = -1
    read(sarg,*) nx
    call error(nx<1, 'nx must be positive integer')

    ! read ny
    call get_command_argument(2, sarg)
    ny = -1
    read(sarg,*) ny
    call error(ny<1, 'ny must be positive integer')

    ! read nt
    call get_command_argument(3, sarg)
    nt = -1
    read(sarg,*) nt
    call error(nt<1, 'nt must be positive integer')

    ! read total time
    call get_command_argument(4, sarg)
    t = -1.
    read(sarg,*) t
    call error(t<0, 't must be positive real value')

    if ( nargs == 5) then
        verbose_output = .true.
    else
        verbose_output = .false.
    end if

    ! store the parameters
    options%global_nx = nx
    options%global_ny = ny
    options%N  = nx*ny
    options%nt = nt
    ! compute timestep size
    options%dt = t/nt
    ! compute the distance between grid points
    ! assume that x dimension has length 1.0
    options%dx = 1./real(nx-1)
    ! set alpha, assume diffusion coefficient D is 1
    options%alpha = (options%dx**2) / (1.*options%dt)
end
!==============================================================================

!initialize MPI
subroutine initialize_mpi(options, domain)
    implicit none

    ! arguments
    type(discretizationT), intent(inout):: options
    type(subdomainT),      intent(out)  :: domain

    ! local
    integer     :: err
    integer     :: mpi_rank, mpi_size
    integer     :: ndomx, ndomy
    integer     :: domx, domy
    integer     :: nx, ny, startx, starty, endx, endy
    integer     :: i

    call mpi_init(err)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)

    ! compute the domain decomposition size
    ! ndomx and ndomy are the number of sub-domains in the x and y directions
    ! repsectively
    ndomx = sqrt(real(mpi_size))
    do while ( mod(mpi_size, ndomx) > 0 )
        ndomx = ndomx-1
    end do
    ndomy = mpi_size / ndomx

    ! compute this sub-domain index
    ! work backwards from: mpi_rank = (domx-1) + (domy-1)*ndomx
    domx = mod(mpi_rank, ndomx) + 1
    domy = (mpi_rank-domx+1)/ndomx + 1

    nx = options%global_nx / ndomx
    ny = options%global_ny / ndomy
    startx = (domx-1)*nx+1
    starty = (domy-1)*ny+1
    ! adjust for grid dimensions that do not divided evenly between the
    ! sub-domains
    if ( domx .eq. ndomx ) then
        nx = options%global_nx - startx + 1
    endif
    if ( domy .eq. ndomy ) then
        ny = options%global_ny - starty + 1
    endif

    endx = startx + nx -1
    endy = starty + ny -1

    domain%startx = startx
    domain%starty = starty
    domain%endx = endx
    domain%endy = endy

    ! store the local grid dimensions back into options
    options%nx = nx
    options%ny = ny
    options%N  = nx*ny

    domain%neighbour_east  = mpi_rank+1
    domain%neighbour_west  = mpi_rank-1
    domain%neighbour_north = mpi_rank+ndomx
    domain%neighbour_south = mpi_rank-ndomx

    if (domx .eq. 1) then
        domain%neighbour_west = -1
    endif
    if (domx .eq. ndomx) then
        domain%neighbour_east = -1
    endif
    if (domy .eq. 1) then
        domain%neighbour_south = -1
    endif
    if (domy .eq. ndomy) then
        domain%neighbour_north = -1
    endif

    domain%rank = mpi_rank
    domain%size = mpi_size

    domain%ndomx = ndomx
    domain%ndomy = ndomy
    domain%domx = domx
    domain%domy = domy

    do i=0,mpi_size-1
        if( mpi_rank == i .and. verbose_output ) then
            write (*,'(A,I2,A,I2,A,I2,A,I2,A,A,I2,A,I2,A,A,I3,I3,A,I3,I3,A,I5,A,I5)') 'rank ', mpi_rank, ' /', mpi_size, &
                        ' : (', domx, ',', domy, ' )', &
                        '/(', ndomx, ',', ndomy, ') ', &
                        ' neigh N-S ', domain%neighbour_north, domain%neighbour_south, &
                        ' neigh E-W ', domain%neighbour_east, domain%neighbour_west, &
                        ' local dims ', nx, '  x' , ny
        end if
        call mpi_barrier(MPI_COMM_WORLD, ierr);
    end do
end

!==============================================================================

! write error message and terminate
subroutine error(yes, msg)
    implicit none

    ! arguments
    logical, intent(in) :: yes
    character(len=*), intent(in) :: msg

    ! local
    integer, external :: lnblnk

    if (yes) then
        write(0,*) 'FATAL PROGRAM ERROR!'
        write(0,*) msg
        write(0,*) 'Execution aborted...'
        stop
    end if
end

!==============================================================================

end program
