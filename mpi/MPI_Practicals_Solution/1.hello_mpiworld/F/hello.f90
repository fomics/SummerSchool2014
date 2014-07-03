PROGRAM hello

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
! Purpose: A simple MPI-program printing "Hello world!"        !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE


  INTEGER ierror


  CALL MPI_INIT(ierror)

  WRITE(*,*) 'Hello world!'

  CALL MPI_FINALIZE(ierror)

END PROGRAM
