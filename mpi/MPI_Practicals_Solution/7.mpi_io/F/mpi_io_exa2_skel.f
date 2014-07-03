        PROGRAM mpi_io_test

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C This file has been written as a sample solution to an exercise in a 
C course given at the HLRS www.hlrs.de . It is made
C freely available with the understanding that every copy of this file
C must include this header and that HLRS takes no responsibility for
C the use of the enclosed teaching material.
C
C Authors:    Rolf Rabenseifner
C
C Contact:    rabenseifner@hlrs.de 
C
C Purpose:    A program to test parallel file I/O with MPI.
C
C Contents:   F source code.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        IMPLICIT NONE

        INCLUDE 'mpif.h'

        INTEGER ierror, my_rank, size, i

        INTEGER ndims, array_of_sizes(1), array_of_subsizes(1)
        INTEGER array_of_starts(1), order
        INTEGER fh
        INTEGER etype
        INTEGER filetype
        INTEGER (____) disp 
        INTEGER status(MPI_STATUS_SIZE) 

        CHARACTER buf 

        CALL MPI_INIT(ierror)

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

        etype = MPI_CHARACTER
        ndims = ____
        array_of_sizes(1)    = ____
        array_of_subsizes(1) = ____
        array_of_starts(1)   = ____
        order = MPI_ORDER_FORTRAN 
        CALL MPI_TYPE_CREATE_SUBARRAY(ndims, array_of_sizes,
     &                               array_of_subsizes, array_of_starts,
     &                               order, etype, filetype, ierror)
        CALL MPI_TYPE____

        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'my_test_file', 
     &                     MPI_MODE_____ + MPI_MODE_____,
     &                     MPI_INFO_NULL, fh, ierror) 
 
        disp = ____
        CALL MPI_FILE_SET_VIEW(fh, disp, etype, filetype, 'native',
     &                         MPI_INFO_NULL, ierror) 

        DO I=1,3
          buf = CHAR( ICHAR('a') + my_rank ) 
          CALL MPI_FILE_WRITE(fh, buf, ____, ____, status, ierror) 
        END DO 

        CALL MPI_FILE_CLOSE(fh, ierror)

        WRITE(*,*) 'PE=',my_rank

        CALL MPI_FINALIZE(ierror)

        STOP
        END
