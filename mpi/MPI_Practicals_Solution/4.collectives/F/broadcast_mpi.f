      program mpi_reduction
      implicit none
      include 'mpif.h'
      integer rank, data, ierror
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      if (rank.eq.0) then
         write (*,*) 'enter value'
         read (*,*) data
      end if
! broadcast the value of data of rank 0 to all ranks
      call MPI_Bcast(data, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
      write (*,*) "I am rank", rank, "and the value is", data
      call MPI_Finalize(ierror)
      end program mpi_reduction
