
!  |-----------|
!  | 0| 4| 8|12|
!  |-----------|
!  | 1| 5| 9|13|
!  |-----------|
!  | 2| 6|10|14|
!  |-----------|
!  | 3| 7|11|15|
!  |-----------|

!  xggggggggggx
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  gddddddddddg
!  xggggggggggx

      program ghost_cell_exchange
      implicit none
      include 'mpif.h'
      integer SUBDOMAIN, DOMAINSIZE
      parameter (SUBDOMAIN = 10)
      parameter (DOMAINSIZE = SUBDOMAIN+2)
      integer rank, size, i, j, rank_right, rank_left, ierror
      integer rank_top, rank_bottom
      integer dims(2), periods(2)
      double precision data(DOMAINSIZE,DOMAINSIZE)
      integer request
      INTEGER status(MPI_STATUS_SIZE)
      integer comm_cart
      integer data_ghost
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
      if (size.ne.16) then
         write (*,*)"please run this with 16 processors"
         call MPI_Finalize(ierror)
         stop
      end if
      do i=1, DOMAINSIZE, 1
         do j=1, DOMAINSIZE, 1
            data(i,j)=rank
         end do
      end do
! neighbouring ranks with cartesian grid communicator
      dims(1)=4
      dims(2)=4
      periods(1)=1
      periods(2)=1
! we do not allow the reordering of ranks here
! an alternative solution would be to allow the reordering and to use the new communicator for the communication
! then the MPI library has the opportunity to choose the best rank order with respect to performance
      call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0,
     &                     comm_cart, ierror)
      call MPI_Cart_shift(comm_cart, 0, 1, rank_left, rank_right,
     &                    ierror)
      call MPI_Cart_shift(comm_cart, 1, 1, rank_top, rank_bottom,
     &                    ierror)
! derived datatype
      call MPI_Type_vector(SUBDOMAIN, 1, DOMAINSIZE, MPI_DOUBLE,
     &                     data_ghost, ierror)
      call MPI_Type_commit(data_ghost, ierror)
!  ghost cell exchange with the neighbouring cells in all directions
!  a) MPI_Send, MPI_Irecv
!  b) MPI_Isend, MPI_Recv
!  c) MPI_Sendrecv
!  to the left
      call MPI_Irecv(data(2,DOMAINSIZE), SUBDOMAIN, MPI_DOUBLE,
     &               rank_right, 0, MPI_COMM_WORLD, request, ierror)
      call MPI_Send(data(2,2), SUBDOMAIN, MPI_DOUBLE, rank_left, 0,
     &              MPI_COMM_WORLD, ierror)
      call MPI_Wait(request, status, ierror)
!  to the right
      call MPI_Irecv(data(2,1), SUBDOMAIN, MPI_DOUBLE, rank_left, 0,
     &               MPI_COMM_WORLD, request, ierror)
      call MPI_Send(data(2,DOMAINSIZE-1), SUBDOMAIN, MPI_DOUBLE,
     &              rank_right, 0, MPI_COMM_WORLD, ierror)
      call MPI_Wait(request, status, ierror)
!  to the top
      call MPI_Irecv(data(DOMAINSIZE,2), 1, data_ghost, rank_bottom, 0,
     &               MPI_COMM_WORLD, request, ierror)
      call MPI_Send(data(2,2), 1, data_ghost, rank_top, 0,
     &              MPI_COMM_WORLD, ierror)
      call MPI_Wait(request, status, ierror)
!  to the bottom
      call MPI_Irecv(data(1,2), 1, data_ghost, rank_top, 0,
     &               MPI_COMM_WORLD, request, ierror)
      call MPI_Send(data(SUBDOMAIN+1,2), 1, data_ghost, rank_bottom, 0,
     &              MPI_COMM_WORLD, ierror)
      call MPI_Wait(request, status, ierror)
      if (rank.eq.4) then
         write (*,*) 'data of rank 4 after communication'
         do i=1, DOMAINSIZE, 1
            do j=1, DOMAINSIZE, 1
               write (*,'(F6.1)',advance='no') data(i,j)
            end do
            write (*,*)
         end do
      end if
      call MPI_Type_free(data_ghost)
      call MPI_Comm_free(comm_cart)
      call MPI_Finalize(ierror)
      end program ghost_cell_exchange
