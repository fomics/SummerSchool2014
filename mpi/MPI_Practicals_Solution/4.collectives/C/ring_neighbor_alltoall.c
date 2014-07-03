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
 * Purpose: Creating a 1-dimensional topology.                  *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define to_right 201
#define max_dims 1 


int main (int argc, char *argv[])
{
  int my_rank, size;
  int snd_buf_arr[2], rcv_buf_arr[2];
   /* snd_buf = snd_buf_arr[1],  rcv_buf = rcv_buf_arr[0] */

  int right, left;
  int sum, i;

  MPI_Comm    new_comm;
  int         dims[max_dims],
              periods[max_dims],
              reorder;
/*int         my_coords[max_dims]; */

  MPI_Status  status;
  MPI_Request request;


  MPI_Init(&argc, &argv);

  /* Get process info. */
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Set cartesian topology. */
  dims[0] = size;
  periods[0] = 1;
  reorder = 1;
 
  MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, periods,
                  reorder,&new_comm);

  /* Get coords */
  MPI_Comm_rank(new_comm, &my_rank);
  /* MPI_Cart_coords(new_comm, my_rank, max_dims, my_coords); */ 

  /* Get nearest neighbour rank. */
  MPI_Cart_shift(new_comm, 0, 1, &left, &right);

  /* Compute global sum. */
 
  sum = 0;
  snd_buf_arr[1] = my_rank;

  for( i = 0; i < size; i++) 
  {
    /* MPI_Issend(&snd_buf_arr[1], 1, MPI_INT, right, to_right, new_comm, &request);
       MPI_Recv(&rcv_buf_arr[0], 1, MPI_INT, left, to_right, new_comm, &status);
       MPI_Wait(&request, &status);
    */
    MPI_Neighbor_alltoall(snd_buf_arr, 1, MPI_INT, rcv_buf_arr, 1, MPI_INT, new_comm);
    
    snd_buf_arr[1] = rcv_buf_arr[0];
    sum += rcv_buf_arr[0];
  }

  printf ("PE%i:\tSum = %i\n", my_rank, sum);
  /* printf ("PE%i, Coords = %i: Sum = %i\n", 
              my_rank, my_coords[0], sum); */

  MPI_Finalize();
}
