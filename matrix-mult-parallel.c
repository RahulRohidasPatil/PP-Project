#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_SIZE 744

void allocateMatrix(int ***, int);
void fillMatrix(int **, int);
void printMatrix(int **, int);
void freeMatrix(int **);
void Local_matrix_allocate(int ***, int, int);
void Local_matrix_multiply(int **, int **, int **, int, int);

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int size, rank, **matrixA, **matrixB, **matrixC, **localMatrixA, **localMatrixC;
    int elements_count = MATRIX_SIZE * MATRIX_SIZE;
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (MATRIX_SIZE % size)
    {
        if (rank == 0)
            printf("Matrix Size must be divisible by number of processors.\n");
        MPI_Finalize();
        return 0;
    }

    allocateMatrix(&matrixA, MATRIX_SIZE);
    allocateMatrix(&matrixB, MATRIX_SIZE);
    allocateMatrix(&matrixC, MATRIX_SIZE);

    if (rank == 0)
    {
        fillMatrix(matrixA, MATRIX_SIZE);
        printf("Matrix A:\n");
        // printMatrix(matrixA, MATRIX_SIZE);

        fillMatrix(matrixB, MATRIX_SIZE);
        printf("Matrix B:\n");
        // printMatrix(matrixB, MATRIX_SIZE);
    }

    MPI_Bcast(*matrixB, elements_count, MPI_INT, 0, MPI_COMM_WORLD);

    int sendCount = elements_count / size;
    int received_row_count = sendCount / MATRIX_SIZE;

    Local_matrix_allocate(&localMatrixA, received_row_count, MATRIX_SIZE);
    MPI_Scatter(*matrixA, sendCount, MPI_INT, *localMatrixA, sendCount, MPI_INT, 0, MPI_COMM_WORLD);

    Local_matrix_allocate(&localMatrixC, received_row_count, MATRIX_SIZE);
    Local_matrix_multiply(localMatrixA, matrixB, localMatrixC, received_row_count, MATRIX_SIZE);

    MPI_Gather(*localMatrixC, sendCount, MPI_INT, *matrixC, sendCount, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("Matrix C:\n");
        // printMatrix(matrixC, MATRIX_SIZE);

        freeMatrix(matrixA);
        freeMatrix(matrixB);
        freeMatrix(matrixC);
        freeMatrix(localMatrixA);
        freeMatrix(localMatrixC);
    }

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}

void allocateMatrix(int ***matrix, int size)
{
    *matrix = (int **)malloc(size * sizeof(int *));
    int *block = (int *)malloc(size * size * sizeof(int));
    for (int i = 0; i < size; i++)
        (*matrix)[i] = block + i * size;
}

void fillMatrix(int **matrix, int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[i][j] = rand() % 10;
}

void printMatrix(int **matrix, int size)
{
    printf("[");
    for (int i = 0; i < size; i++)
    {
        printf("%s", i == 0 ? "[" : " [");
        for (int j = 0; j < size; j++)
            printf("%d%s", matrix[i][j], j != size - 1 ? ", " : "");
        printf("%s", i != size - 1 ? "]\n" : "]");
    }
    printf("]\n");
}

void freeMatrix(int **matrix)
{
    free(matrix[0]);
    free(matrix);
}

void Local_matrix_allocate(int ***matrix, int rows, int columns)
{
    *matrix = (int **)malloc(rows * sizeof(int *));
    int *block = (int *)malloc(rows * columns * sizeof(int));
    for (int i = 0; i < rows; i++)
        (*matrix)[i] = block + i * columns;
}

void Local_matrix_multiply(int **matrixA, int **matrixB, int **matrixC, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
        {
            matrixC[i][j] = 0;
            for (int k = 0; k < columns; k++)
                matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
        }
}
