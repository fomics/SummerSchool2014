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
 * Purpose: A program to meassure 1-dim halo communication      *
 *          in myrank -1 and +1 directions (left and right)     *
 *          using MPI_Neighbor_alltoall                         *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

#include <stdio.h>
#include <mpi.h>

#define number_of_messages 50
#define start_length 4
#define length_factor 8
#define max_length 8388608 /* ==> 2 x 32 MB per process */
#define number_package_sizes 8
/* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
/* #define number_package_sizes 9 */
#define max_dims 1 

int main(int argc, char *argv[])
{
  int i, j, length, my_rank, left, right, size, test_value, mid;    
  double start, finish, transfer_time; 
/*  MPI_Request rq[2]; */
/*  MPI_Status status_arr[2]; */
  float snd_buf[2*max_length];
  float rcv_buf[2*max_length];
/*      snd_buf_left[i] = snd_buf[i] ;  snd_buf_right[i] = snd_buf[i+length] */
/*      rcv_buf[i] = rcv_buf[i] ;  rcv_buf[i+length] = rcv_buf[i+length] */

  MPI_Comm    new_comm;
  int         dims[max_dims],
              periods[max_dims],
              reorder;

  /* Get process info. */
/* Naming conventions                                                                */
/* Processes:                                                                        */
/*     my_rank-1                        my_rank                         my_rank+1    */
/* "left neighbor"                     "myself"                     "right neighbor" */
/*   ...    rcv_buf_right <--- snd_buf_left snd_buf_right ---> rcv_buf_left    ...   */
/*   ... snd_buf_right ---> rcv_buf_left       rcv_buf_right <--- snd_buf_left ...   */
/*                        |                                  |                       */
/*              halo-communication                 halo-communication                */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Set cartesian topology. */
  dims[0] = size;
  periods[0] = 1;
  reorder = 1;
 
  MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, periods,
                  reorder,&new_comm);

  /* Get my_rank gain in reordered new_comm */
  MPI_Comm_rank(new_comm, &my_rank);
  /* Get nearest neighbour rank. */
  MPI_Cart_shift(new_comm, 0, 1, &left, &right);

  if (my_rank == 0) printf("    message size      transfertime  duplex bandwidth per process and neighbor\n");

  length = start_length;

  for (j = 1; j <= number_package_sizes; j++)
  { 
    for (i = 0; i <= number_of_messages; i++)
    {
      if(i==1) start = MPI_Wtime();

      test_value = j*1000000 + i*10000 + my_rank*10 ; mid = (length-1)/number_of_messages*i;

      snd_buf[0]=test_value+1  ; snd_buf[mid]=test_value+2  ; snd_buf[length-1]=test_value+3;
      snd_buf[0+length]=test_value+6 ; snd_buf[mid+length]=test_value+7 ; snd_buf[length-1+length]=test_value+8;

/*      MPI_Irecv(rcv_buf+length, length, MPI_FLOAT, right, 17, new_comm, &rq[0]); */
/*      MPI_Irecv(rcv_buf,        length, MPI_FLOAT, left,  23, new_comm, &rq[1]); */
/*      MPI_Send(snd_buf,        length, MPI_FLOAT, left,  17, new_comm); */
/*      MPI_Send(snd_buf+length, length, MPI_FLOAT, right, 23, new_comm); */
/*      MPI_Waitall(2, rq, status_arr); */

      MPI_Neighbor_alltoall(snd_buf, length, MPI_FLOAT, rcv_buf, length, MPI_FLOAT, new_comm); 

/*    ...snd_buf_... is used to store the values that were stored in snd_buf_... in the neighbor process */
      test_value = j*1000000 + i*10000 + left*10  ; mid = (length-1)/number_of_messages*i;
      snd_buf[0+length]=test_value+6 ; snd_buf[mid+length]=test_value+7 ; snd_buf[length-1+length]=test_value+8;
      test_value = j*1000000 + i*10000 + right*10 ; mid = (length-1)/number_of_messages*i;
      snd_buf[0]=test_value+1  ; snd_buf[mid]=test_value+2  ; snd_buf[length-1]=test_value+3;
      if ((rcv_buf[0] != snd_buf[0+length]) || (rcv_buf[mid] != snd_buf[mid+length]) || 
                                                   (rcv_buf[length-1] != snd_buf[length-1+length])) {
         printf("%d: j=%d, i=%d --> snd_buf[0,%d,%d+length]=(%f,%f,%f)\n",
                    my_rank, j, i, mid, length-1, snd_buf[0+length], snd_buf[mid+length], snd_buf[length-1+length]);
         printf("%d:     is not identical to rcv_buf[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank,       mid, length-1, rcv_buf[0],  rcv_buf[mid],  rcv_buf[length-1]);
      }
      if ((rcv_buf[0+length] != snd_buf[0]) || (rcv_buf[mid+length] != snd_buf[mid]) ||
                                                   (rcv_buf[length-1+length] != snd_buf[length-1])) {
         printf("%d: j=%d, i=%d --> snd_buf[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank, j, i, mid, length-1, snd_buf[0],  snd_buf[mid],  snd_buf[length-1]);
         printf("%d:     is not identical to rcv_buf[%d,%d,%d+length]=(%f,%f,%f)\n",
                    my_rank,    0+length, mid+length, length-1+length, rcv_buf[0+length], rcv_buf[mid+length], rcv_buf[length-1+length]);
      }

    }
    finish = MPI_Wtime();

    if (my_rank == 0) 
    {
      transfer_time = (finish - start) / number_of_messages;
      printf("%10i bytes %12.3f usec %13.3f MB/s\n", 
             length*(int)sizeof(float), transfer_time*1e6, 1.0e-6*2*length*sizeof(float) / transfer_time);
    }

    length = length * length_factor;
  }

  MPI_Finalize();
}
