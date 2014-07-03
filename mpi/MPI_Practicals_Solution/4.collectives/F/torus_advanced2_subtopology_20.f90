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
! Purpose: Creating a 2-dimensional topology.                  !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: max_dims=2

  INTEGER :: ierror, my_rank, size

  INTEGER :: sum

  INTEGER :: new_comm, slice_comm
  INTEGER :: dims(max_dims)
  LOGICAL :: reorder, periods(max_dims), remain_dims(max_dims)
  INTEGER :: coords(max_dims), size_of_slice, my_rank_in_slice 


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

! Set two-dimensional cartesian topology.
  dims(1) = 0       
  dims(2) = 0
  periods(1) = .TRUE.
  periods(2) = .TRUE.
  reorder = .TRUE.

  CALL MPI_DIMS_CREATE(size, max_dims, dims, ierror)
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, max_dims, dims, &
                           periods, reorder, new_comm, ierror)
  CALL MPI_COMM_RANK(new_comm, my_rank, ierror)
  CALL MPI_CART_COORDS(new_comm,my_rank, max_dims, coords, ierror) 

! Split the new-comm into slices
  remain_dims(1) = .TRUE.
  remain_dims(2) = .FALSE.
  CALL MPI_CART_SUB(new_comm, remain_dims, slice_comm, ierror)
  CALL MPI_COMM_SIZE(slice_comm, size_of_slice, ierror)
  CALL MPI_COMM_RANK(slice_comm, my_rank_in_slice, ierror)

! Compute sum.

  CALL MPI_ALLREDUCE(my_rank, sum, 2, MPI_INTEGER, MPI_SUM, slice_comm, ierror)

  WRITE(*,*) "PE",my_rank, ", coords = (",coords(1),",", coords(2), &
                 "), Slice_rank=", my_rank_in_slice, ": Sum =", sum

  CALL MPI_FINALIZE(ierror)

END PROGRAM
