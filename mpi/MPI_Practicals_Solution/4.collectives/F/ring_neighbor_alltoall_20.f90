PROGRAM ring

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
! Purpose: Using MPI_Neighbor_alltoall for ring communication. !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201

  INTEGER, PARAMETER :: max_dims=1

  INTEGER :: ierror, my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER :: snd_buf_arr(0:1)
!            snd_buf = snd_buf_arr(1)
  INTEGER :: rcv_buf_arr(0:1)
!            rcv_buf = rcv_buf_arr(0)

  INTEGER :: status(MPI_STATUS_SIZE)

  INTEGER :: request

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy

  INTEGER :: new_comm          
  INTEGER :: dims(max_dims)
  LOGICAL :: reorder, periods(max_dims)
! INTEGER :: coords(max_dims) 


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

! Set one-dimensional cartesian topology.

  dims(1) = size
  periods(1) = .TRUE.
  reorder = .TRUE.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD, max_dims, dims, &
                           periods, reorder, new_comm, ierror)
  CALL MPI_COMM_RANK(new_comm, my_rank, ierror)
! CALL MPI_CART_COORDS(new_comm,my_rank, max_dims, coords, ierror) 

! Get nearest neighbour ranks.

  CALL MPI_CART_SHIFT(new_comm, 0, 1, left, right, ierror)

! Compute sum.

  sum = 0
  snd_buf_arr(1) = my_rank

  DO i = 1, size

!      CALL MPI_ISSEND(snd_buf_arr(1), 1, MPI_INTEGER, right, to_right, new_comm, request, ierror)
!      CALL MPI_RECV(rcv_buf_arr(0), 1, MPI_INTEGER, left, to_right, new_comm, status, ierror)
!      CALL MPI_WAIT(request, status, ierror)

     CALL MPI_NEIGHBOR_ALLTOALL(snd_buf_arr, 1, MPI_INTEGER, rcv_buf_arr, 1, MPI_INTEGER, new_comm, ierror)

!      CALL MPI_GET_ADDRESS(snd_buf_arr, iadummy, ierror)
! !    ... or with MPI-3.0 and later:
! !    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf_arr)

     snd_buf_arr(1) = rcv_buf_arr(0)
     sum = sum + rcv_buf_arr(0)

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum =", sum
! WRITE(*,*) "PE", my_rank, " Coord =", coords(1), ": Sum =", sum

  CALL MPI_FINALIZE(ierror)

END PROGRAM
