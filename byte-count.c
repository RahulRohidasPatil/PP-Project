#include <stdio.h>
#include <string.h>
#include <mpi.h>

const int MAX_STRING = 100;

int main(int argc, char *argv[])
{
    char greeting[MAX_STRING];
    int size, rank, count;
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    if (rank != 0)
    {
        sprintf(greeting, "Greetings from process %d of %d!", rank, size);
        MPI_Send(greeting, strlen(greeting) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        printf("Greetings from process %d of %d!\n", rank, size);
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(greeting, MAX_STRING, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("%s\n", greeting);
            MPI_Get_count(&status, MPI_CHAR, &count);
            printf("byte-count: %d\n", count);
        }
    }

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}