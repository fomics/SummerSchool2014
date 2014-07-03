PROGRAM ring_derived

!==============================================================!
!                                                              !
! This file has been written as a sample solution to an        !
! exercise in a course given at the High Performance           !
! Computing Centre Stuttgart (HLRS).                           !
! The examples are based on the examples in the MPI course of  !
! the Edinburgh Parallel Computing Centre (EPCC).              !
! It is made freely available with the understanding that      !
! every copy of this file must include this header and that    !
! HLRS and EPCC take no responsibility for the use of the      !
! enclosed teaching material.                                  !
!                                                              !
! Authors: Joel Malard, Alan Simpson,            (EPCC)        !
!          Rolf Rabenseifner, Traugott Streicher (HLRS)        !
!                                                              !
! Contact: rabenseifner@hlrs.de                                !
!                                                              !
! Purpose: A program that uses derived data-types.             !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201

  INTEGER :: ierror, my_rank, size

  INTEGER :: right, left

  INTEGER :: i

  TYPE t
     SEQUENCE
     INTEGER :: i
     REAL    :: r
  END TYPE t
  TYPE(t) :: sum, snd_buf, rcv_buf

  INTEGER(KIND=MPI_ADDRESS_KIND) :: first_var_address, second_var_address
  INTEGER :: send_recv_type

  INTEGER :: array_of_block_length(2)
  INTEGER :: array_of_types(2)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_displacements(2)

  INTEGER :: status(MPI_STATUS_SIZE)


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

  right = mod(my_rank+1,      size)
  left  = mod(my_rank-1+size, size)
!     ... this SPMD-style neighbor computation with modulo has the same meaning as:
!     right = my_rank + 1
!     IF (right .EQ. size) right = 0
!     left = my_rank - 1
!     IF (left .EQ. -1) left = size-1

! Create derived datatype. 

  array_of_block_length(1) = 1
  array_of_block_length(2) = 1

  array_of_types(1) = MPI_INTEGER
  array_of_types(2) = MPI_REAL

  CALL MPI_GET_ADDRESS(snd_buf%i, first_var_address, ierror)
  CALL MPI_GET_ADDRESS(snd_buf%r, second_var_address, ierror)

  array_of_displacements(1) = 0
  array_of_displacements(2) = second_var_address - first_var_address

  CALL MPI_TYPE_CREATE_STRUCT(2, array_of_block_length, array_of_displacements, array_of_types, send_recv_type, ierror)
  CALL MPI_TYPE_COMMIT(send_recv_type, ierror)

  sum%i = 0
  sum%r = 0
  snd_buf%i = my_rank
  snd_buf%r = REAL(my_rank)

  DO i = 1, size

     CALL MPI_SENDRECV(snd_buf, 1, send_recv_type, right, to_right,  &
                       rcv_buf, 1, send_recv_type, left, to_right,   &
                       MPI_COMM_WORLD, status, ierror)

     snd_buf = rcv_buf
     sum%i = sum%i + rcv_buf%i
     sum%r = sum%r + rcv_buf%r

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum%i =", sum%i, " Sum%r =", sum%r

  CALL MPI_FINALIZE(ierror)

END PROGRAM
