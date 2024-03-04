#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_SIZE 256

void allocateMatrix(int ***, int);
void fillMatrix(int **, int);
void printMatrix(int **, int);
void freeMatrix(int **);
void multiplyMatrices(int **, int **, int **, int);

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int size, rank, **matrixA, **matrixB, **matrixC;
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        allocateMatrix(&matrixA, MATRIX_SIZE);
        fillMatrix(matrixA, MATRIX_SIZE);
        printf("Matrix A:\n");
        // printMatrix(matrixA, MATRIX_SIZE);

        allocateMatrix(&matrixB, MATRIX_SIZE);
        fillMatrix(matrixB, MATRIX_SIZE);
        printf("Matrix B:\n");
        // printMatrix(matrixB, MATRIX_SIZE);

        allocateMatrix(&matrixC, MATRIX_SIZE);

        multiplyMatrices(matrixA, matrixB, matrixC, MATRIX_SIZE);
        printf("Matrix C:\n");
        // printMatrix(matrixC, MATRIX_SIZE);
    }

    if (rank == 0)
    {
        freeMatrix(matrixA);
        freeMatrix(matrixB);
        freeMatrix(matrixC);
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

void multiplyMatrices(int **matrixA, int **matrixB, int **matrixC, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrixC[i][j] = 0;
            for (int k = 0; k < size; k++)
                matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
        }
    }
}
