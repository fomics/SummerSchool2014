#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]){
   int rank, size, number;
   char string_comm[1000];
   MPI_Status stat;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   if (size!=2){
      printf("please run this with 2 processors\n");
      MPI_Finalize();
      exit(1);
   }
   if (rank==0){
      printf("enter number \n");
      scanf("%d",&number);
   }
// send the contents of number from rank 0 to rank 1 using MPI_Send --- MPI_Recv




   if (rank==1)
      printf("number communicated %i\n", number);
   MPI_Finalize();
   return(0);
}
