        program exa_block 
        implicit none 

        include 'mpif.h'

        integer i, ierror
        integer gsize, lsize, psize, rank, nblocks
        parameter (gsize=10000000)
 
!       real*8 garray(gsize)
!HPF$   PROCESSORS   P(psize)
!HPF$   DISTRIBUTE garray(BLOCK)
!       real*8 larray(lsize)

!       inline function definition for 
!       division with rounding to next upper integer 
        integer   ii,jj, UPPER_DIV 
        UPPER_DIV(ii,jj) = (ii+jj-1) / jj

        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, psize, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank,  ierror) 
 
        lsize = UPPER_DIV(gsize, psize) 
        nblocks = UPPER_DIV(gsize, lsize) 
        if(rank.GE.nblocks)   lsize = 0 
        if(rank.EQ.nblocks-1) lsize = gsize - (nblocks-1)*lsize
         
!       write(*,*)'rank=',rank,'  g=',gsize,'  l=',lsize,'  p=',psize
        if (rank.eq.0) write(*,*) 'gsize=',gsize, ' on',psize,'PEs'

        call sub(rank,gsize,lsize,psize)

        end 

!-----------------------------------------------------------------------

        subroutine sub(rank,gsize,lsize,psize)
        implicit none 
        include 'mpif.h'
        integer rank,gsize,lsize,psize
 
        integer i, darray_type, ierror, fh
        integer distribs, dargs, istatus(MPI_STATUS_SIZE)
        integer (kind=MPI_OFFSET_KIND) disp 
 
        real*8 larray(lsize)
        integer lext 

        double precision start_time, end_time 

        distribs = MPI_DISTRIBUTE_BLOCK
        dargs    = MPI_DISTRIBUTE_DFLT_DARG

        call MPI_TYPE_CREATE_DARRAY(psize,rank,1,gsize,distribs,
     &                               dargs, psize, MPI_ORDER_FORTRAN,
     &                               MPI_REAL8, darray_type,ierror)
        call MPI_TYPE_COMMIT(darray_type,ierror) 

        call MPI_TYPE_EXTENT(MPI_REAL8, lext, ierror) 
 
        do i=1,lsize
         larray(i) = 100*i + rank
        enddo
 
!       --------------- 
 
        call MPI_BARRIER(MPI_COMM_WORLD, ierror) ! only for benchmarking
        start_time = MPI_WTIME() 

        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'exa_block_testfile',
     &                     MPI_MODE_CREATE + MPI_MODE_WRONLY, 
     &                     MPI_INFO_NULL, fh, ierror)
        disp = 0 
        call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, darray_type,
     &                     'native'    ,MPI_INFO_NULL,ierror)
 
        call MPI_FILE_WRITE_ALL(fh, larray, lsize, MPI_REAL8, 
     &                                       istatus,ierror)

        call MPI_FILE_CLOSE(fh,ierror)
 
        call MPI_BARRIER(MPI_COMM_WORLD, ierror) ! only for benchmarking
        end_time = MPI_WTIME() 

!       write(*,*)'rank=',rank,'  done'
 
        call MPI_BARRIER(MPI_COMM_WORLD, ierror) ! only for benchmarking
        if (rank.eq.0) then
          write(*,*) 'WRITE_ALL in', end_time - start_time, 'sec',
     &        ' ==> ', 1e-6*gsize*lext/(end_time-start_time), 'MB/s'  
        endif 

!       --------------- 
 
        call MPI_FINALIZE(ierror)
 
        return 
        end 
