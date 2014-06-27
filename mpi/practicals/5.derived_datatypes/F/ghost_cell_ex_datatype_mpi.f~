
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
      double precision data(DOMAINSIZE,DOMAINSIZE)
      integer request
      INTEGER status(MPI_STATUS_SIZE)
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
      rank_right=mod(rank+4,16)
      rank_left=mod(rank+16-4,16)
      rank_top=mod(rank+16-1,4)+(rank/4)*4
      rank_bottom=mod(rank+1,4)+(rank/4)*4
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
      call MPI_Finalize(ierror)
      end program ghost_cell_exchange
