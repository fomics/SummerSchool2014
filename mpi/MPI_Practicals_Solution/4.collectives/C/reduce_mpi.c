#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
   int rank, input, result;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   printf("enter value\n");
//   scanf ("%d",&input);
   input=rank+1;
// reduce the values of the different ranks in input to result of rank 0 with the operation sum (max, logical and)
   MPI_Reduce(&input, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
// reduce the values of the different ranks in input to result of all ranks with the operation sum (max, logical and)
   MPI_Allreduce(&input, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank==0)
      printf("result %i\n", result);
   MPI_Finalize();
   return(0);
}
