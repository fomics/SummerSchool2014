PROGRAM halo

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
! Purpose: A program to meassure 1-dim halo communication      !
!          in myrank -1 and +1 directions (left and right)     !
!          using MPI_Neighbor_alltoall                         !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi
  IMPLICIT NONE

  INTEGER, PARAMETER :: number_of_messages=50
  INTEGER, PARAMETER :: start_length=4
  INTEGER, PARAMETER :: length_factor=8
  INTEGER, PARAMETER :: max_length=8388608     ! ==> 2 x 32 MB per process
  INTEGER, PARAMETER :: number_package_sizes=8
! INTEGER, PARAMETER :: max_length=67108864    ! ==> 2 x 0.5 GB per process
! INTEGER, PARAMETER :: number_package_sizes=9

  INTEGER, PARAMETER :: max_dims=1

  INTEGER :: i, j, length, my_rank, left, right, size, test_value, mid, ierror
  INTEGER(KIND=MPI_ADDRESS_KIND) :: lb, size_of_real
  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy
  DOUBLE PRECISION :: start, finish, transfer_time
! INTEGER :: rq(2) 
! INTEGER :: status_arr(MPI_STATUS_SIZE,2)
! REAL :: snd_buf_left(max_length), snd_buf_right(max_length)
! REAL, ASYNCHRONOUS :: rcv_buf_left(max_length), rcv_buf_right(max_length)
  REAL :: snd_buf(max_length*2)
! REAL, ASYNCHRONOUS :: rcv_buf(max_length*2)
  REAL :: rcv_buf(max_length*2)
!         snd_buf_left(i) = snd_buf(i) ;  snd_buf_right(i) = snd_buf(i+length)
!         rcv_buf_left(i) = rcv_buf(i) ;  rcv_buf_right(i) = rcv_buf(i+length)

  INTEGER :: new_comm          
  INTEGER :: dims(max_dims)
  LOGICAL :: reorder, periods(max_dims)

! Naming conventions
! Processes: 
!     my_rank-1                        my_rank                         my_rank+1
! "left neighbor"                     "myself"                     "right neighbor"
!   ...    rcv_buf_right <--- snd_buf_left snd_buf_right ---> rcv_buf_left    ...
!   ... snd_buf_right ---> rcv_buf_left       rcv_buf_right <--- snd_buf_left ...
!                        |                                  |
!              halo-communication                 halo-communication

  CALL MPI_INIT(ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)

! Set one-dimensional cartesian topology.

  dims(1) = size
  periods(1) = .TRUE.
  reorder = .TRUE.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD, max_dims, dims, &
                           periods, reorder, new_comm, ierror)
  CALL MPI_COMM_RANK(new_comm, my_rank, ierror)

! Get nearest neighbour ranks.

  CALL MPI_CART_SHIFT(new_comm, 0, 1, left, right, ierror)

  CALL MPI_TYPE_GET_EXTENT(MPI_REAL, lb, size_of_real, ierror) 

  IF (my_rank .EQ. 0) THEN
     WRITE (*,*) "message size   transfertime    duplex bandwidth per process and neigbor"
  END IF

  length = start_length

  DO j = 1, number_package_sizes

     DO i = 0, number_of_messages

        IF (i==1) start = MPI_WTIME()

        test_value = j*1000000 + i*10000 + my_rank*10 ; mid = 1 + (length-1)/number_of_messages*i

        snd_buf(1)=test_value+1 ; snd_buf(mid)=test_value+2 ; snd_buf(length)=test_value+3
        snd_buf(1+length)=test_value+6 ; snd_buf(mid+length)=test_value+7 ; snd_buf(length+length)=test_value+8

!         CALL MPI_IRECV(rcv_buf(1+length), length, MPI_REAL, right, 17, new_comm, rq(1), ierror)
!         CALL MPI_IRECV(rcv_buf(1),        length, MPI_REAL, left,  23, new_comm, rq(2), ierror)
!         CALL MPI_SEND (snd_buf(1),        length, MPI_REAL, left,  17, new_comm, ierror)
!         CALL MPI_SEND (snd_buf(1+length), length, MPI_REAL, right, 23, new_comm, ierror)
!         CALL MPI_WAITALL(2, rq, status_arr, ierror)
!         CALL MPI_GET_ADDRESS(rcv_buf,  iadummy, ierror)
! !       ... or with MPI-3.0 and later:
! !       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)

        CALL MPI_NEIGHBOR_ALLTOALL(snd_buf, length, MPI_REAL, rcv_buf, length, MPI_REAL, new_comm, ierror)

!       ...snd_buf_... is used to store the values that were stored in snd_buf_... in the neighbor process
        test_value = j*1000000 + i*10000 + left*10  ; mid = 1 + (length-1)/number_of_messages*i
        snd_buf(1+length)=test_value+6 ; snd_buf(mid+length)=test_value+7 ; snd_buf(length+length)=test_value+8
        test_value = j*1000000 + i*10000 + right*10 ; mid = 1 + (length-1)/number_of_messages*i
        snd_buf(1)=test_value+1 ; snd_buf(mid)=test_value+2 ; snd_buf(length)=test_value+3
        IF ((rcv_buf(1).NE.snd_buf(1+length)).OR.(rcv_buf(mid).NE.snd_buf(mid+length)).OR. &
                                                     (rcv_buf(length).NE.snd_buf(length+length))) THEN
           write (*,*) my_rank,": j=",j," i=",i," -->  snd_buf(",1,mid,length,"+length)=", &
                                                     snd_buf(1+length),snd_buf(mid+length),snd_buf(length+length)," in left process"
           write (*,*) my_rank,":     is not identical to rcv_buf(",1,mid,length,")=", &
                                                     rcv_buf(1),rcv_buf(mid),rcv_buf(length) 
        END IF
        IF ((rcv_buf(1+length).NE.snd_buf(1)).OR.(rcv_buf(mid+length).NE.snd_buf(mid)).OR. &
                                                     (rcv_buf(length+length).NE.snd_buf(length))) THEN
           write (*,*) my_rank,": j=",j," i=",i," <-- snd_buf(",1,mid,length,")=", &
                                                     snd_buf(1),snd_buf(mid),snd_buf(length)," in right process"
           write (*,*) my_rank,":     is not identical to rcv_buf(",1,mid,length,"+length)=", &
                                                     rcv_buf(1+length),rcv_buf(mid+length),rcv_buf(length+length) 
        END IF

     END DO
     finish = MPI_WTIME()
     
     IF (my_rank .EQ. 0) THEN
        transfer_time = (finish - start) / (number_of_messages)
        WRITE(*,*) INT(length*size_of_real),'bytes  ', transfer_time*1e6,'usec  ', 1e-6*2*length*size_of_real/transfer_time,'MB/s'
     END IF

     length = length * length_factor
  END DO

  CALL MPI_FINALIZE(ierror)
END PROGRAM
