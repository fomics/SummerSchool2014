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
 * Purpose: Benchmarking MPI_Ssend and MPI_Recv.                *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define process_A 0
#define process_B 1
#define ping  17
#define pong  23

#define number_of_messages 50
#define start_length 8
#define length_faktor 64
#define max_length 2097152     /* 2 Mega */
#define number_package_sizes 4

int main(int argc, char *argv[])
{
  int my_rank;
 
  int i, j;  
  int length_of_message;    
  double start, finish, time, transfer_time; 
  MPI_Status status;
  float buffer[max_length];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == process_A) 
  {
    printf("message size\ttransfertime\t\tbandwidth\n");
  }
  length_of_message = start_length;

  for (i = 1; i <= number_package_sizes; i++)
  { 
    
    if (my_rank == process_A) 
    {
        MPI_Send(buffer, length_of_message, MPI_FLOAT, process_B, 
                          ping, MPI_COMM_WORLD);

        MPI_Recv(buffer, length_of_message, MPI_FLOAT, process_B, 
                       pong, MPI_COMM_WORLD, &status);
    }
    else if (my_rank == process_B) 
    {
        MPI_Recv(buffer, length_of_message, MPI_FLOAT, process_A, 
                         ping, MPI_COMM_WORLD, &status);
      
        MPI_Send(buffer, length_of_message, MPI_FLOAT, process_A, 
                          pong, MPI_COMM_WORLD);
    }  

    start = MPI_Wtime();
    for (j = 1; j <= number_of_messages; j++)
    {
      if (my_rank == process_A) 
      {
        MPI_Send(buffer, length_of_message, MPI_FLOAT, process_B, 
                          ping, MPI_COMM_WORLD);

        MPI_Recv(buffer, length_of_message, MPI_FLOAT, process_B, 
                       pong, MPI_COMM_WORLD, &status);
      }
      else if (my_rank == process_B) 
      {
        MPI_Recv(buffer, length_of_message, MPI_FLOAT, process_A, 
                         ping, MPI_COMM_WORLD, &status);
      
        MPI_Send(buffer, length_of_message, MPI_FLOAT, process_A, 
                          pong, MPI_COMM_WORLD);
      }  
    }
    finish = MPI_Wtime();

    if (my_rank == process_A) 
    {
      time = finish - start;

      transfer_time = time / (2 * number_of_messages);

      printf("%i bytes\t\t%f usec\t\t%f MB/s\n", 
             length_of_message*(int)sizeof(float),
             transfer_time*1e6,
             1.0e-6*length_of_message*sizeof(float) / transfer_time);
    }
    length_of_message *= length_faktor;

  }
  MPI_Finalize();
}
