#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
   int rank, data;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if (rank==0){
      printf("enter value\n");
      scanf ("%d",&data);
   }
// broadcast the value of data of rank 0 to all ranks
   MPI_Bcast(&data, 1, MPI_INT, 0, MPI_COMM_WORLD);
   printf("I am rank %i and the value is %i\n", rank, data);
   MPI_Finalize();
   return(0);
}
