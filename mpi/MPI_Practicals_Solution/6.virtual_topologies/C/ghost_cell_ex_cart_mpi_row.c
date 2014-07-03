//  C style two-dimensional array rowwise
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//|-----------|
//| 0| 1| 2| 3|
//|-----------|
//| 4| 5| 6| 7|
//|-----------|
//| 8| 9|10|11|
//|-----------|
//|12|13|14|15|
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
   int rank, size, i, j, dims[2], periods[2], rank_top, rank_bottom, rank_left, rank_right;
   double data[DOMAINSIZE*DOMAINSIZE];
   MPI_Request request;
   MPI_Status status;
   MPI_Comm comm_cart;
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
   dims[0]=dims[1]=4;
   periods[0]=periods[1]=1;
//  we do not allow the reordering of ranks here
//  an alternative solution would be to allow the reordering and to use the new communicator for the communication
//  then the MPI library has the opportunity to choose the best rank order with respect to performance
   MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);
   MPI_Cart_shift(comm_cart, 0, 1, &rank_top, &rank_bottom);
   MPI_Cart_shift(comm_cart, 1, 1, &rank_left, &rank_right);
//  derived datatype
   MPI_Type_vector(SUBDOMAIN, 1, DOMAINSIZE, MPI_DOUBLE, &data_ghost);
   MPI_Type_commit(&data_ghost);
//  ghost cell exchange with the neighbouring cells in all directions
//  a) MPI_Send, MPI_Irecv
//  b) MPI_Isend, MPI_Recv
//  c) MPI_Sendrecv
//  to the top
   MPI_Irecv(&data[2-1+(DOMAINSIZE-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_bottom, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[2-1+(2-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_top, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
//  to the bottom
   MPI_Irecv(&data[2-1+(1-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_top, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[2-1+(DOMAINSIZE-1-1)*DOMAINSIZE], SUBDOMAIN, MPI_DOUBLE, rank_bottom, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
//  to the left
   MPI_Irecv(&data[DOMAINSIZE-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_right, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[2-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_left, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
//  to the right
   MPI_Irecv(&data[1-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_left, 0, MPI_COMM_WORLD, &request);
   MPI_Send(&data[DOMAINSIZE-1-1+(2-1)*DOMAINSIZE], 1, data_ghost, rank_right, 0, MPI_COMM_WORLD);
   MPI_Wait(&request, &status);
   if (rank==4){
      printf("data of rank 4 after communication\n");
      for (j=0; j<DOMAINSIZE; j++){
         for (i=0; i<DOMAINSIZE; i++){
            printf("%f ", data[i+j*DOMAINSIZE]);
         }
         printf("\n");
      }
   }
   MPI_Type_free(&data_ghost);
   MPI_Comm_free(&comm_cart);
   MPI_Finalize();
   return(0);
}
