!// process decomposition on 4*4 grid 
!  |-----------|
!  | 0| 4| 8|12|
!  |-----------|
!  | 1| 5| 9|13|
!  |-----------|
!  | 2| 6|10|14|
!  |-----------|
!  | 3| 7|11|15|
!  |-----------|

!Each process works on a 10*10 (SUBDOMAIN) block of data                                                                               
! the d corresponds to data, g corresponds to "ghost cells"                                                                            
! and x are empty (not exchanged for now)      

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


! task: each rank has to find its left and right neighbor and                                                                          
! send them the data they need (left array goes to left neighbor                                                                       ! and right array goes to bottom neighbor)
 
      program ghost_cell_exchange
      implicit none
      include 'mpif.h'
      integer SUBDOMAIN, DOMAINSIZE
      parameter (SUBDOMAIN = 10)
      parameter (DOMAINSIZE = SUBDOMAIN+2)
      integer rank, size, i, j, rank_right, rank_left, ierror
      double precision data(DOMAINSIZE,DOMAINSIZE)
      integer request
      INTEGER status(MPI_STATUS_SIZE)
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
      if (size.ne.16) then
         write (*,*)"please run this with 16 processors"
         call MPI_Finalize(ierror);
         stop
      end if
      do i=1, DOMAINSIZE, 1
         do j=1, DOMAINSIZE, 1
            data(i,j)=rank
         end do
      end do
      rank_right=!find the rank of my right neighbour
      rank_left=!find the rank of my left neighbour

!  ghost cell exchange with the neighbouring cells on the left and on the right
!  a) MPI_Send, MPI_Irecv
!  b) MPI_Isend, MPI_Recv
!  c) MPI_Sendrecv
!  to the left
! a)

! b)

! c)

!  to the right
! a)

! b)

! c)

         write (*,*) 'data of rank 4 after communication'
         do i=1, DOMAINSIZE, 1
            do j=1, DOMAINSIZE, 1
               write (*,'(F6.1)',advance='no') data(i,j)
            end do
            write (*,*)
         end do
      end if
      call MPI_Finalize(ierror);
      end program ghost_cell_exchange
