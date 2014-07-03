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
 * Purpose: A program with derived datatypes.                   *
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

  struct buff{
     int   i;
     float f;
  } snd_buf, rcv_buf, sum;

  int i;

  int          array_of_blocklengths[COUNT];
  MPI_Aint     array_of_displacements[COUNT], first_var_address, second_var_address;
  MPI_Datatype array_of_types[COUNT], datatype;

  MPI_Status  status;
  MPI_Request request;


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

  MPI_Get_address(&snd_buf.i, &first_var_address);
  MPI_Get_address(&snd_buf.f, &second_var_address);

  array_of_displacements[0] = (MPI_Aint) 0;
  array_of_displacements[1] = second_var_address - first_var_address;

  array_of_types[0] = MPI_INT;
  array_of_types[1] = MPI_FLOAT;

  MPI_Type_create_struct(COUNT, array_of_blocklengths, array_of_displacements, array_of_types, &datatype);
  MPI_Type_commit(&datatype);

  /* Compute global sum. */
  sum.i = 0;            sum.f = 0;
  snd_buf.i = my_rank;  snd_buf.f = my_rank;  /* Step 1 = init */

  for( i = 0; i < size; i++) 
  {
    MPI_Issend(&snd_buf, 1, datatype, right, to_right, MPI_COMM_WORLD, &request);  /* Step 2a */
    MPI_Recv(&rcv_buf, 1, datatype, left, to_right, MPI_COMM_WORLD, &status);      /* Step 3 */
    MPI_Wait(&request, &status);                                                  /* Step 2b */
    snd_buf = rcv_buf;                        /* Step 4 */
    sum.i += rcv_buf.i;  sum.f += rcv_buf.f;  /* Step 5 */
  }

  printf ("PE%i:\tSum = %i\t%f\n", my_rank, sum.i, sum.f);

  MPI_Finalize();
}
