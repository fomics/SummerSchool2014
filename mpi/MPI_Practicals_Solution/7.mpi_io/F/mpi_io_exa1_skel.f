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

        INTEGER fh
        INTEGER (____) offset
        INTEGER status(MPI_STATUS_SIZE) 

        CHARACTER buf 

        CALL MPI_INIT(ierror)

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'my_test_file', 
     &                     MPI_MODE_____ + MPI_MODE_____,
     &                     MPI_INFO_NULL, fh, ierror) 
 
        DO I=1,10
          buf = CHAR( ICHAR('0') + my_rank ) 
          offset = ________________ 
          CALL MPI_FILE_WRITE_AT(fh, offset, buf, ____, ____,
     &                           status, ierror) 
        END DO 

        CALL MPI_FILE_CLOSE(fh, ierror)

        WRITE(*,*) 'PE=',my_rank

        CALL MPI_FINALIZE(ierror)

        STOP
        END
