//  FORTRAN style two-dimensional array columnwise
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//|-----------|
//| 0| 4| 8|12|
//|-----------|
//| 1| 5| 9|13|
//|-----------|
//| 2| 6|10|14|
//|-----------|
//| 3| 7|11|15|
//|-----------|

// xggggggggggx
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// gddddddddddg
// xggggggggggx

#define SUBDOMAIN 10
#define DOMAINSIZE (SUBDOMAIN+2)

int main(int argc, char *argv[]){
   int rank, size, i, j, rank_left, rank_right, rank_top, rank_bottom;
   double data[DOMAINSIZE*DOMAINSIZE];
   MPI_Request request;
   MPI_Status status;
   MPI_Datatype data_ghost;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   if (size!=16){
      printf("please run this with 16 processors\n");
      MPI_Finalize();
      exit(1);
   }
   for (i=0; i<DOMAINSIZE*DOMAINSIZE; i++){
      data[i]=rank;
   }
//  neighbouring ranks with cartesian grid communicator
   rank_left=(rank+16-4)%16;
   rank_right=(rank+4)%16;
   rank_top=(rank+16-1)%4+(rank/4)*4;
   rank_bottom=(rank+1)%4+(rank/4)*4;
//  derived datatype
   MPI_Type_vector(SUBDOMAIN, 1, DOMAINSIZE, MPI_DOUBLE, &data_ghost);
   MPI_Type_commit(&data_ghost);
//  ghost cell exchange with the neighbouring cells in all directions
//  a) MPI_Send, MPI_Irecv
//  b) MPI_Isend, MPI_Recv
//  c) MPI_Sendrecv
//  to the left
   MPI_Irecv(&data[2-1+(DOMAINSIZE-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[2-1+(2-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
//  to the right
   MPI_Irecv(&data[2-1+(1-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_left, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[2-1+(DOMAINSIZE-1-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_right, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
//  to the top
   MPI_Irecv(&data[DOMAINSIZE-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_bottom, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[2-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_top, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
//  to the bottom
   MPI_Irecv(&data[1-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_top, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[DOMAINSIZE-1-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_bottom, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
   if (rank==4){
      printf("data of rank 4 after communication\n");
      for (i=0; i<DOMAINSIZE; i++){
         for (j=0; j<DOMAINSIZE; j++){
            printf("%f ", data[i+j*DOMAINSIZE]);
         }
         printf("\n");
      }
   }
   MPI_Type_free(&data_ghost);
   MPI_Finalize();
   return(0);
}
