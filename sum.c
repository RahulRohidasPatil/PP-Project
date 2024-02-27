#include <stdio.h>
#include <string.h>
#include <mpi.h>

const int MAX_STRING = 100;

int main(int argc, char *argv[])
{
    int size, rank, sum = 0;
    int numbers[48] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0)
    {
        int startIndex = rank * 2;
        int local_sum = numbers[startIndex] + numbers[startIndex + 1];
        MPI_Send(&local_sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        sum += numbers[0] + numbers[1];
        for (int i = 1; i < size; i++)
        {
            int recv_sum;
            MPI_Recv(&recv_sum, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += recv_sum;
        }
        printf("Sum: %d\n", sum);
    }

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}