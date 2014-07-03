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
 * Purpose: A program to try MPI_Comm_size and MPI_Comm_rank.   *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

#include <stdio.h>
#include <mpi.h>


int main(int argc, char *argv[])
{
   int my_rank, size;


   MPI_Init(&argc, &argv);

   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (my_rank == 0)   
   {
     printf ("Hello world!\n");
   }

   printf("I am process %i out of %i.\n", my_rank, size);

   MPI_Finalize();
}
