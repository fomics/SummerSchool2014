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
 * Purpose: A program to try MPI_Ssend and MPI_Recv.            *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define proc_A 0
#define proc_B 1
#define ping  17
#define pong  23
#define number_of_messages 50
#define length_of_message   1

int main(int argc, char *argv[])
{
  int my_rank;

  float buffer[length_of_message];

  int i;   

  double start, finish, time;

  MPI_Status status;


  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  start = MPI_Wtime();
  for (i = 1; i <= number_of_messages; i++)
  {
    if (my_rank == proc_A) 
    {
      MPI_Send(buffer, length_of_message, MPI_FLOAT, proc_B, ping,
                       MPI_COMM_WORLD);
      MPI_Recv(buffer, length_of_message, MPI_FLOAT, proc_B, pong,
                       MPI_COMM_WORLD, &status);
    }
    else if (my_rank == proc_B) 
    { 
      MPI_Recv(buffer, length_of_message, MPI_FLOAT, proc_A, ping,
                       MPI_COMM_WORLD, &status);
      MPI_Send(buffer, length_of_message, MPI_FLOAT, proc_A, pong,
                       MPI_COMM_WORLD);
    }
  }
  finish = MPI_Wtime();

  if (my_rank == proc_A) 
  {
    time = finish - start;
    printf("Time for one messsage: %f micro seconds.\n", 
           (float)(time / (2 * number_of_messages) * 1e6));
  } 

  MPI_Finalize();
}
