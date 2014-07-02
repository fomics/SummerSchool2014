module io

! dependencies
use mpi
use data,   only: subdomainT, discretizationT, options, domain

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_header(filename)
    ! arguments
    character(len=*)          :: filename

    ! variables
    integer :: output
    output=20

    ! metadata
    open (unit=output, file='output.bov', status='replace')
    write(output,*) 'TIME: 0.0'
    write(output,*) 'DATA_FILE: output.bin'
    write(output,*) 'DATA_SIZE: ', options%global_nx, ' ', options%global_ny, ' 1'
    write(output,*) 'DATA_FORMAT: DOUBLE'
    write(output,*) 'VARIABLE: phi'
    write(output,*) 'DATA_ENDIAN: LITTLE'
    write(output,*) 'CENTERING: nodal'
    write(output,*) 'BYTE_OFFSET: 0' ! using MPI IO the byte offset is 0
    write(output,*) 'BRICK_SIZE: ', 1.0  , ' ', real(options%global_ny-1)*options%dx , ' 1.0'
    close (output)
end subroutine write_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_parallel(filename, u)
    ! arguments
    real (kind=8), intent(in) :: u(options%nx,options%ny)
    character(len=*)          :: filename

    ! variables
    integer(kind=mpi_offset_kind) :: disp
    integer :: ierr, filehandle, filetype, N
    integer :: dimuids(2), ustart(2), ucount(2)

    ! initial displacement is zero
    disp = 0

    ! open file handle
    call mpi_file_open(                         &
        MPI_COMM_WORLD, Filename,               &
        ior(MPI_MODE_CREATE, MPI_MODE_WRONLY),  &
        MPI_INFO_NULL, filehandle, ierr)

    ! calculate field dimensions
    ustart(1) = domain%startx-1
    ustart(2) = domain%starty-1

    ucount(1) = options%nx
    ucount(2) = options%ny
    N = ucount(1) * ucount(2)

    dimuids(1) = options%global_nx
    dimuids(2) = options%global_ny

    ! write header
    !if (domain%rank == 0) then
        !call mpi_file_write(filehandle, dimuids, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    !endif
    !disp = sizeof(dimuids)

    ! create a subarray representing the local block
    call MPI_Type_create_subarray(2, dimuids, ucount, ustart, &
                                       MPI_ORDER_FORTRAN, MPI_DOUBLE, filetype, ierr)
    call MPI_Type_commit(filetype, ierr)

    call MPI_File_set_view(filehandle, disp, MPI_DOUBLE, &
                           filetype, "native", MPI_INFO_NULL, ierr)

    call MPI_File_write_all(filehandle, u, N, MPI_DOUBLE, MPI_STATUS_IGNORE, ierr)
    call MPI_Type_free(filetype, ierr)
    call MPI_File_close(filehandle, ierr)
end subroutine write_parallel

end module io
