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
 * Purpose: A program to try MPI_Sendrecv_replace.              *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define to_right 201
#define COUNT 2


int main (int argc, char *argv[])
{
  int my_rank, size;
  int right, left;

  int   int_buf, int_sum, i;
  float float_buf, float_sum;

  int          array_of_blocklengths[COUNT];
  MPI_Aint     array_of_displacements[COUNT], first_var_address, second_var_address;
  MPI_Datatype array_of_types[COUNT], sendrecvtype;

  MPI_Status  status;

  /* Get process and neighbour info. */
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  right = (my_rank+1)      % size;
  left  = (my_rank-1+size) % size;
/* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
/* right = my_rank + 1;          */
/* if (right == size) right = 0; */
/* left = my_rank - 1;           */
/* if (left == -1) left = size-1;*/

  /* Set MPI datatypes for sending and receiving partial sums. */

  array_of_blocklengths[0] = 1;
  array_of_blocklengths[1] = 1;

  MPI_Address(&int_buf, &first_var_address);
  MPI_Address(&float_buf, &second_var_address);

  array_of_displacements[0] = (MPI_Aint) 0;
  array_of_displacements[1] = second_var_address - first_var_address;

  array_of_types[0] = MPI_INT;
  array_of_types[1] = MPI_FLOAT;

  MPI_Type_struct(COUNT, array_of_blocklengths, array_of_displacements, array_of_types, &sendrecvtype);

  MPI_Type_commit(&sendrecvtype);

  /* Compute global sum. */
  int_sum = 0;
  float_sum = 0;
  int_buf = my_rank;
  float_buf = (float) my_rank;

  for( i = 0; i < size; i++) 
  {
    MPI_Sendrecv_replace(&int_buf, 1, sendrecvtype, right, to_right, 
                                                    left,  to_right, MPI_COMM_WORLD, &status) ;
    
    int_sum += int_buf;
    float_sum += float_buf;
  }

  printf ("PE%i:\tSum = %i\t%f\n", my_rank, int_sum, float_sum);

  MPI_Finalize();
}
