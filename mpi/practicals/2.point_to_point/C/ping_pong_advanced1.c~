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
 *modified by Themis Athanassiadou                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define proc_A 0
#define proc_B 1
#define ping  17
#define pong  23
#define number_of_messages 50
#define length_of_message   1

/* Modify your ping-pong code such that you exchange one message before entering the 50 
   message loop*/

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

  finish = MPI_Wtime();

  if (my_rank == proc_A) 
  {
    time = finish - start;
    printf("Time for one messsage: %f micro seconds.\n", 
           (float)(time / (2 * number_of_messages) * 1e6));
  } 

  MPI_Finalize();
}
