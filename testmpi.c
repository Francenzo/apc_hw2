#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int myid, numprocs;
    int tag, source, destination, count;
    double *sendbuff = (double *)malloc(sizeof(double) * 1);
    double *recvbuff = (double *)malloc(sizeof(double) * 1);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Request request1, request2;
    MPI_Status status;
    /* Write your code here */
    count = 1;
    tag = 0;
    //only two threads
    if (myid == 0)
    {
        //send request
        *sendbuff = 3.14159;
        source = myid;
        destination = 1;

        MPI_Isend(sendbuff, 1, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD, &request1);
        MPI_Wait(&request1, &status);
    }
    else
    {
        //recieve request
        source = 0;
        MPI_Irecv(recvbuff, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &request2);
        MPI_Wait(&request2, &status);
        printf("thread (%d) get the buffer %f\n", myid, *recvbuff);
    }

    MPI_Finalize();
}
