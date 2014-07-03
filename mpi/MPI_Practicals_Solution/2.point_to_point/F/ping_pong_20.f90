PROGRAM pingpong

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
! Purpose: A program to try MPI_Ssend and MPI_Recv.            !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER proc_a
  PARAMETER(proc_a=0)
            
  INTEGER proc_b
  PARAMETER(proc_b=1)                

  INTEGER ping
  PARAMETER(ping=17)
        
  INTEGER pong
  PARAMETER(pong=23)        

  INTEGER number_of_messages 
  PARAMETER (number_of_messages=50)

  INTEGER length
  PARAMETER (length=1)
 
  DOUBLE PRECISION start, finish, time
  INTEGER status(MPI_STATUS_SIZE)
   
  REAL buffer(length)

  INTEGER i

  INTEGER ierror, my_rank, size


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)

  start = MPI_WTIME()

  DO i = 1, number_of_messages

     IF (my_rank .EQ. proc_a) THEN
        CALL MPI_SEND(buffer, length, MPI_REAL, proc_b, ping, MPI_COMM_WORLD, ierror)
        CALL MPI_RECV(buffer, length, MPI_REAL, proc_b, pong, MPI_COMM_WORLD, status, ierror)
     ELSE IF (my_rank .EQ. proc_b) THEN
        CALL MPI_RECV(buffer, length, MPI_REAL, proc_a, ping, MPI_COMM_WORLD, status, ierror)
        CALL MPI_SEND(buffer, length, MPI_REAL, proc_a, pong, MPI_COMM_WORLD, ierror)
     END IF

  END DO

  finish = MPI_WTIME()

  IF (my_rank .EQ. proc_a) THEN

     time = finish - start

     WRITE(*,*) 'Time for one message:', time/(2*number_of_messages)*1e6, ' micro seconds'

  END IF

  CALL MPI_FINALIZE(ierror)

END PROGRAM
