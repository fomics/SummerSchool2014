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
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi
  USE, INTRINSIC ::ISO_C_BINDING
  IMPLICIT NONE

  INTEGER, PARAMETER :: number_of_messages=50
  INTEGER, PARAMETER :: start_length=4
  INTEGER, PARAMETER :: length_factor=8
  INTEGER, PARAMETER :: max_length=8388608     ! ==> 2 x 32 MB per process
  INTEGER, PARAMETER :: number_package_sizes=8
! INTEGER, PARAMETER :: max_length=67108864    ! ==> 2 x 0.5 GB per process
! INTEGER, PARAMETER :: number_package_sizes=9

  INTEGER :: i, j, length, my_rank, left, right, size, test_value, mid, ierror
  INTEGER(KIND=MPI_ADDRESS_KIND) :: lb, size_of_real
  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy
  DOUBLE PRECISION :: start, finish, transfer_time
  REAL :: snd_buf_left(max_length), snd_buf_right(max_length)
  REAL, POINTER, ASYNCHRONOUS :: rcv_buf_left(:), rcv_buf_right(:)

  INTEGER :: win_rcv_buf_left, win_rcv_buf_right
  INTEGER :: disp_unit
  INTEGER(KIND=MPI_ADDRESS_KIND) :: buf_size, target_disp
  TYPE(C_PTR) :: ptr_rcv_buf_left, ptr_rcv_buf_right
  INTEGER :: offset_left, offset_right

  INTEGER :: rq(2) 
  INTEGER :: status_arr(MPI_STATUS_SIZE,2)
  INTEGER :: snd_dummy_left, snd_dummy_right
  INTEGER, ASYNCHRONOUS :: rcv_dummy_left, rcv_dummy_right
  INTEGER :: k

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
  left  = mod(my_rank-1+size, size)
  right = mod(my_rank+1,      size)

  CALL MPI_TYPE_GET_EXTENT(MPI_REAL, lb, size_of_real, ierror) 

  buf_size = max_length * size_of_real
  disp_unit = size_of_real
  CALL MPI_WIN_ALLOCATE_SHARED(buf_size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, ptr_rcv_buf_left, win_rcv_buf_left, ierror)
  CALL C_F_POINTER(ptr_rcv_buf_left, rcv_buf_left, (/max_length/))
! offset_left  is defined so that rcv_buf_left(xxx+offset_left) in process 'my_rank' is the same location as
!                                 rcv_buf_left(xxx) in process 'left':
  offset_left  = +(left-my_rank)*max_length

  CALL MPI_WIN_ALLOCATE_SHARED(buf_size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, ptr_rcv_buf_right, win_rcv_buf_right, ierror)
  CALL C_F_POINTER(ptr_rcv_buf_right, rcv_buf_right, (/max_length/))
! offset_right is defined so that rcv_buf_right(xxx+offset_right) in process 'my_rank' is the same location as
!                                 rcv_buf_right(xxx) in process 'right':
  offset_right  = +(right-my_rank)*max_length

  target_disp = 0

  IF (my_rank .EQ. 0) THEN
     WRITE (*,*) "message size   transfertime    duplex bandwidth per process and neigbor"
  END IF

  length = start_length

  CALL MPI_WIN_LOCK_ALL(MPI_MODE_NOCHECK, win_rcv_buf_left,  ierror)
  CALL MPI_WIN_LOCK_ALL(MPI_MODE_NOCHECK, win_rcv_buf_right, ierror)

  DO j = 1, number_package_sizes

     DO i = 0, number_of_messages

        IF (i==1) start = MPI_WTIME()

        test_value = j*1000000 + i*10000 + my_rank*10 ; mid = 1 + (length-1)/number_of_messages*i

        snd_buf_left(1)=test_value+1  ; snd_buf_left(mid)=test_value+2  ; snd_buf_left(length)=test_value+3
        snd_buf_right(1)=test_value+6 ; snd_buf_right(mid)=test_value+7 ; snd_buf_right(length)=test_value+8

        CALL MPI_GET_ADDRESS(rcv_buf_left,  iadummy, ierror)
        CALL MPI_GET_ADDRESS(rcv_buf_right, iadummy, ierror)
!       ... or with MPI-3.0 and later:
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_left)
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_right)

!       ... The local WIN_SYNCs are not needed to sync the processor and real memory on the target-process after the processors sync
!           because there is no "data transfer" on the shared memory
!       CALL MPI_WIN_SYNC(win_rcv_buf_right, ierror)
!       CALL MPI_WIN_SYNC(win_rcv_buf_left,  ierror)
 
!       ... tag=17: posting to left that rcv_buf_left can be stored from left // tag=23: and rcv_buf_right from right
        CALL MPI_IRECV(rcv_dummy_left,  1, MPI_INTEGER, right, 17, MPI_COMM_WORLD, rq(1), ierror)
        CALL MPI_IRECV(rcv_dummy_right, 1, MPI_INTEGER, left,  23, MPI_COMM_WORLD, rq(2), ierror)
        CALL MPI_SEND (snd_dummy_right, 0, MPI_INTEGER, left,  17, MPI_COMM_WORLD, ierror)
        CALL MPI_SEND (snd_dummy_left,  0, MPI_INTEGER, right, 23, MPI_COMM_WORLD, ierror)
        CALL MPI_WAITALL(2, rq, status_arr, ierror)

!       ... The local WIN_SYNCs are not needed to sync the processor and real memory on the target-process after the processors sync
!           because there is no "data transfer" on the shared memory
!       CALL MPI_WIN_SYNC(win_rcv_buf_right, ierror)
!       CALL MPI_WIN_SYNC(win_rcv_buf_left,  ierror)

        CALL MPI_GET_ADDRESS(rcv_buf_left,  iadummy, ierror)
        CALL MPI_GET_ADDRESS(rcv_buf_right, iadummy, ierror)
!       ... or with MPI-3.0 and later:
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_left)
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_right)

!       CALL MPI_PUT(snd_buf_left,  length, MPI_REAL, left,  target_disp, length, MPI_REAL, win_rcv_buf_right, ierror)
!       CALL MPI_PUT(snd_buf_right, length, MPI_REAL, right, target_disp, length, MPI_REAL, win_rcv_buf_left,  ierror)
!       ... is substited by:
!       rcv_buf_right(1+offset_left  : length+offset_left)  = snd_buf_left(1:length)
!       rcv_buf_left (1+offset_right : length+offset_right) = snd_buf_right(1:length)
        IF (length>1000000) THEN
          DO k = 1, length
            rcv_buf_right(k+offset_left ) = snd_buf_left(k)
          END DO
          DO k = 1, length
            rcv_buf_left (k+offset_right) = snd_buf_right(k)
          END DO
        ELSE
          DO k = 1, length
            rcv_buf_right(k+offset_left ) = snd_buf_left(k)
            rcv_buf_left (k+offset_right) = snd_buf_right(k)
          END DO
        ENDIF

        CALL MPI_GET_ADDRESS(rcv_buf_left,  iadummy, ierror)
        CALL MPI_GET_ADDRESS(rcv_buf_right, iadummy, ierror)
!       ... or with MPI-3.0 and later:
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_left)
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_right)

!       ... The local WIN_SYNCs are needed to sync the processor and real memory on the origin-process before the processors sync
        CALL MPI_WIN_SYNC(win_rcv_buf_right, ierror)
        CALL MPI_WIN_SYNC(win_rcv_buf_left,  ierror)
 
!       ... The following communication synchronizes the processors in the way that the origin processor has finished the store
!            before the target processor starts to load the data.  
!       ... tag=17: posting to right that rcv_buf_left was stored from left // tag=23: and rcv_buf_right from right
        CALL MPI_IRECV(rcv_dummy_left,  1, MPI_INTEGER, left,  17, MPI_COMM_WORLD, rq(1), ierror)
        CALL MPI_IRECV(rcv_dummy_right, 1, MPI_INTEGER, right, 23, MPI_COMM_WORLD, rq(2), ierror)
        CALL MPI_SEND (snd_dummy_right, 0, MPI_INTEGER, right, 17, MPI_COMM_WORLD, ierror)
        CALL MPI_SEND (snd_dummy_left,  0, MPI_INTEGER, left,  23, MPI_COMM_WORLD, ierror)
        CALL MPI_WAITALL(2, rq, status_arr, ierror)

!       ... The local WIN_SYNCs are needed to sync the processor and real memory on the target-process after the processors sync
        CALL MPI_WIN_SYNC(win_rcv_buf_right, ierror)
        CALL MPI_WIN_SYNC(win_rcv_buf_left,  ierror)

        CALL MPI_GET_ADDRESS(rcv_buf_left,  iadummy, ierror)
        CALL MPI_GET_ADDRESS(rcv_buf_right, iadummy, ierror)
!       ... or with MPI-3.0 and later:
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_left)
!       IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf_right)

!       ...snd_buf_... is used to store the values that were stored in snd_buf_... in the neighbor process
        test_value = j*1000000 + i*10000 + left*10  ; mid = 1 + (length-1)/number_of_messages*i
        snd_buf_right(1)=test_value+6 ; snd_buf_right(mid)=test_value+7 ; snd_buf_right(length)=test_value+8
        test_value = j*1000000 + i*10000 + right*10 ; mid = 1 + (length-1)/number_of_messages*i
        snd_buf_left(1)=test_value+1  ; snd_buf_left(mid)=test_value+2  ; snd_buf_left(length)=test_value+3
        IF ((rcv_buf_left(1).NE.snd_buf_right(1)).OR.(rcv_buf_left(mid).NE.snd_buf_right(mid)).OR. &
                                                     (rcv_buf_left(length).NE.snd_buf_right(length))) THEN
           write (*,*) my_rank,": j=",j," i=",i," -->  snd_buf_right(",1,mid,length,")=", &
                                                     snd_buf_right(1),snd_buf_right(mid),snd_buf_right(length)," in left process"
           write (*,*) my_rank,":     is not identical to rcv_buf_left(",1,mid,length,")=", &
                                                     rcv_buf_left(1),rcv_buf_left(mid),rcv_buf_left(length) 
        END IF
        IF ((rcv_buf_right(1).NE.snd_buf_left(1)).OR.(rcv_buf_right(mid).NE.snd_buf_left(mid)).OR. &
                                                     (rcv_buf_right(length).NE.snd_buf_left(length))) THEN
           write (*,*) my_rank,": j=",j," i=",i," <-- snd_buf_left(",1,mid,length,")=", &
                                                     snd_buf_left(1),snd_buf_left(mid),snd_buf_left(length)," in right process"
           write (*,*) my_rank,":     is not identical to rcv_buf_right(",1,mid,length,")=", &
                                                     rcv_buf_right(1),rcv_buf_right(mid),rcv_buf_right(length) 
        END IF

     END DO
     finish = MPI_WTIME()
     
     IF (my_rank .EQ. 0) THEN
        transfer_time = (finish - start) / (number_of_messages)
        WRITE(*,*) INT(length*size_of_real),'bytes  ', transfer_time*1e6,'usec  ', 1e-6*2*length*size_of_real/transfer_time,'MB/s'
     END IF

     length = length * length_factor
  END DO

  CALL MPI_WIN_UNLOCK_ALL(win_rcv_buf_left,  ierror)
  CALL MPI_WIN_UNLOCK_ALL(win_rcv_buf_right, ierror)

  CALL MPI_WIN_FREE(win_rcv_buf_left, ierror)
  CALL MPI_WIN_FREE(win_rcv_buf_right, ierror)

  CALL MPI_FINALIZE(ierror)
END PROGRAM
