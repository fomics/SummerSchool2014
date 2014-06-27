      program mpi_reduction
      implicit none
      include 'mpif.h'
      integer i, rank, size, senddata(20), receivedata, ierror
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
      if (size.gt.20) then
         if (rank.eq.0) then
            write (*,*) "do not use more than 20 processors"
            call MPI_Finalize(ierror)
         end if
      end if
      if (rank.eq.0) then
         do i=1, size, 1
            write (*,*) 'enter value'
            read (*,*) senddata(i)
         end do
      end if
! scatter the value of senddata of rank 0 to receivedata of all ranks
      call MPI_Scatter(senddata, 1, MPI_INT, receivedata, 1, MPI_INT,
     &                 0, MPI_COMM_WORLD, ierror)
      write (*,*) "I am rank", rank, "and the value is", receivedata
      call MPI_Finalize(ierror)
      end program mpi_reduction
