/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the High Performance           *
 * Computing Centre Stuttgart (HLRS).                           *
 * The examples are based on the examples in the MPI course of  *
 * the Edinburgh Parallel Computing Centre (EPCC).              *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * HLRS and EPCC take no responsibility for the use of the      *
 * enclosed teaching material.                                  *
 *                                                              *
 * Authors: Joel Malard, Alan Simpson,            (EPCC)        *
 *          Rolf Rabenseifner, Traugott Streicher (HLRS)        *
 *                                                              *
 * Contact: rabenseifner@hlrs.de                                * 
 *                                                              *  
 * Purpose: Trying MPI_Scan in a ring topology.                 *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define tag_ready 7781


int main (int argc, char *argv[])
{
  int my_rank, size;
  int sum;
  int token;

  MPI_Status status;


  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Compute partial rank sum. */

  MPI_Scan (&my_rank, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD ); 

  /* Output in natural order */
  /* CAUTION: Although the printing is initialized by the 
              MPI processes in the order of the ranks,
              it is not guaranteed that the merge of the stdout
              of the processes will keep this order
              on the common output of all MPI processes ! */

  if (my_rank != 0)   
  {
    MPI_Recv(&token, 1, MPI_INT, my_rank - 1, tag_ready, MPI_COMM_WORLD, &status);
  }

  printf ("PE%d:\tSum = %d\n", my_rank, sum);  
 
  if (my_rank != size - 1)
  {
    MPI_Send(&token, 1, MPI_INT, my_rank + 1, tag_ready, MPI_COMM_WORLD);
  }

  MPI_Finalize();
}
