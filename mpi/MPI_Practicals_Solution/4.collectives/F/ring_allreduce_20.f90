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
! Purpose: A program to try MPI_Allreduce.                     !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER :: ierror, my_rank

  INTEGER :: sum   

  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)

  CALL MPI_ALLREDUCE(my_rank, sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)

  WRITE(*,*) "PE", my_rank, ": Sum =", sum

  CALL MPI_FINALIZE(ierror)

END PROGRAM
