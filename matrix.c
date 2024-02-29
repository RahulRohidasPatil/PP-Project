#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_SIZE 2

void allocateMatrix(int ***matrix);
void fillMatrix(int **matrix);
void printMatrix(int **matrix);
void freeMatrix(int **matrix);

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int size, rank, **matrixA, **matrixB;
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        allocateMatrix(&matrixA);
        fillMatrix(matrixA);
        printMatrix(matrixA);

        allocateMatrix(&matrixB);
        fillMatrix(matrixB);
        printMatrix(matrixB);
    }

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}

void allocateMatrix(int ***matrix)
{
    *matrix = (int **)malloc(MATRIX_SIZE * sizeof(int *));
    for (int i = 0; i < MATRIX_SIZE; i++)
        (*matrix)[i] = (int *)malloc(MATRIX_SIZE * sizeof(int));
}

void fillMatrix(int **matrix)
{
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++)
            matrix[i][j] = rand() % 10;
}

void printMatrix(int **matrix)
{
    printf("[");
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        printf("%s", i == 0 ? "[" : " [");
        for (int j = 0; j < MATRIX_SIZE; j++)
            printf("%d%s", matrix[i][j], j != MATRIX_SIZE - 1 ? ", " : "");
        printf("%s", i != MATRIX_SIZE - 1 ? "]\n" : "]");
    }
    printf("]\n");
}

void freeMatrix(int **matrix)
{
    for (int i = 0; i < MATRIX_SIZE; i++)
        free(matrix[i]);
    free(matrix);
}