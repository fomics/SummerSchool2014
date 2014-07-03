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
 * Purpose: Creating a 2-dimensional topology.                  *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define to_right 201
#define max_dims 2


int main (int argc, char *argv[])
{
  int my_rank, size;
  int sum;

  MPI_Comm    new_comm, slice_comm;
  int         dims[max_dims],
              periods[max_dims],
              reorder,
              my_coords[max_dims],
              remain_dims[max_dims],
              size_of_slice, my_rank_in_slice;


  MPI_Init(&argc, &argv);

  /* Get process info. */
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Set cartesian topology. */
  dims[0] = 0;       dims[1] = 0;
  periods[0] = 1;    periods[1] = 1;
  reorder = 1;
 
  MPI_Dims_create(size, max_dims, dims);
  MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, periods,
                  reorder,&new_comm);

  /* Get coords */
  MPI_Comm_rank(new_comm, &my_rank);
  MPI_Cart_coords(new_comm, my_rank, max_dims, my_coords); 

  /* Split the new_comm into slices */
  remain_dims[0]=1;  remain_dims[1]=0;
  MPI_Cart_sub(new_comm, remain_dims, &slice_comm);
  MPI_Comm_size(slice_comm, &size_of_slice);
  MPI_Comm_rank(slice_comm, &my_rank_in_slice);

  /* Compute "global" sums for each slice */
  MPI_Allreduce (&my_rank, &sum, 1, MPI_INT, MPI_SUM, slice_comm);

  printf ("PE%i, Coords=(%i,%i), Slice_rank=%i: Sum = %i\n", 
          my_rank, my_coords[0], my_coords[1], my_rank_in_slice, sum);

  MPI_Finalize();
}
