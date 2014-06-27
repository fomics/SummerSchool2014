#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]){
   int i, rank, size, senddata[20], receivedata;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   if (size>20){
      if (rank==0)
         printf("do not use more than 20 processors\n");
      MPI_Finalize();
      exit(1);
   }
   if (rank==0){
      for (i=0; i<size; i++){
         printf("enter value\n");
         scanf ("%d",&senddata[i]);
      }
   }
// scatter the value of senddata of rank 0 to receivedata of all ranks
   MPI_Scatter(&senddata, 1, MPI_INT, &receivedata, 1, MPI_INT, 0, MPI_COMM_WORLD);
   printf("I am rank %i and the value is %i\n", rank, receivedata);
   MPI_Finalize();
   return(0);
}
