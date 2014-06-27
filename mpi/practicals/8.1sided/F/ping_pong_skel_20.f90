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
! Purpose: A program to meassure MPI_Send and MPI_Recv.        !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi
  IMPLICIT NONE

  INTEGER, PARAMETER :: number_of_messages=50
  INTEGER, PARAMETER :: start_length=1
  INTEGER, PARAMETER :: length_factor=8
! INTEGER, PARAMETER :: max_length=134217728
! INTEGER, PARAMETER :: number_package_sizes=10
  INTEGER, PARAMETER :: max_length=16777216
  INTEGER, PARAMETER :: number_package_sizes=9

  INTEGER :: i, j, length, my_rank, size, test_value, mid, ierror
  INTEGER(KIND=MPI_ADDRESS_KIND) :: lb, size_of_real
  DOUBLE PRECISION :: start, finish, time, transfer_time
  INTEGER :: status(MPI_STATUS_SIZE)
  REAL :: snd_buf(max_length), rcv_buf(max_length)

  CALL MPI_INIT(ierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_TYPE_GET_EXTENT(MPI_REAL, lb, size_of_real, ierror) 

  IF (my_rank .EQ. 0) THEN
     WRITE (*,*) "message size   transfertime    bandwidth"
  END IF

  length = start_length

  DO j = 1, number_package_sizes

     DO i = 0, number_of_messages

        IF(i==1) start = MPI_WTIME()

        test_value = j*1000000 + i*1000 ; mid = 1 + (length-1)/number_of_messages*i
        IF (my_rank .EQ. 0) THEN

           snd_buf(1)=test_value+1 ; snd_buf(mid)=test_value+2 ; snd_buf(length)=test_value+3

           CALL MPI_SEND(snd_buf, length, MPI_REAL, 1, 17, MPI_COMM_WORLD, ierror)

           CALL MPI_RECV(rcv_buf, length, MPI_REAL, 1, 23, MPI_COMM_WORLD, status, ierror)

           snd_buf(1)=test_value+5 ; snd_buf(mid)=test_value+6 ; snd_buf(length)=test_value+7
           IF ((rcv_buf(1).NE.snd_buf(1)).OR.(rcv_buf(mid).NE.snd_buf(mid)).OR.(rcv_buf(length).NE.snd_buf(length))) THEN
              write (*,*) "PE 0: error: j=",j," i=",i," pong snd_buf(",1,mid,length,")=",snd_buf(1),snd_buf(mid),snd_buf(length)
              write (*,*) "PE 0:     ... is not identical to rcv_buf(",1,mid,length,")=",rcv_buf(1),rcv_buf(mid),rcv_buf(length) 
           END IF

        ELSE IF (my_rank .EQ. 1) THEN

           CALL MPI_RECV(rcv_buf, length, MPI_REAL, 0, 17, MPI_COMM_WORLD, status, ierror)

           snd_buf(1)=test_value+1 ; snd_buf(mid)=test_value+2 ; snd_buf(length)=test_value+3
           IF ((rcv_buf(1).NE.snd_buf(1)).OR.(rcv_buf(mid).NE.snd_buf(mid)).OR.(rcv_buf(length).NE.snd_buf(length))) THEN
              write (*,*) "PE 1: error: j=",j," i=",i," ping snd_buf(",1,mid,length,")=",snd_buf(1),snd_buf(mid),snd_buf(length)
              write (*,*) "PE 1:     ... is not identical to rcv_buf(",1,mid,length,")=",rcv_buf(1),rcv_buf(mid),rcv_buf(length) 
           END IF

           snd_buf(1)=test_value+5 ; snd_buf(mid)=test_value+6 ; snd_buf(length)=test_value+7

           CALL MPI_SEND(snd_buf, length, MPI_REAL, 0, 23, MPI_COMM_WORLD, ierror)

        END IF
     END DO
     finish = MPI_WTIME()
     
     IF (my_rank .EQ. 0) THEN
        time = finish - start
        transfer_time = time / (2 * number_of_messages)
        WRITE(*,*) INT(length*size_of_real),'bytes  ', transfer_time*1e6,'usec  ', 1e-6*length*size_of_real/transfer_time,'MB/s'
     END IF

     length = length * length_factor
  END DO

  CALL MPI_FINALIZE(ierror)
END PROGRAM
