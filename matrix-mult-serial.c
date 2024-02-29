#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_SIZE 576

void fillMatrices(int matrix1[][MATRIX_SIZE], int matrix2[][MATRIX_SIZE]);
void printMatrix(int matrix[][MATRIX_SIZE]);
void multiplyMatrices(int matrix1[][MATRIX_SIZE], int matrix2[][MATRIX_SIZE], int productMatrix[][MATRIX_SIZE]);

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int size, rank;
    int matrix1[MATRIX_SIZE][MATRIX_SIZE];
    int matrix2[MATRIX_SIZE][MATRIX_SIZE];
    int productMatrix[MATRIX_SIZE][MATRIX_SIZE];
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        fillMatrices(matrix1, matrix2);

        // printMatrix(matrix1);
        // printMatrix(matrix2);

        multiplyMatrices(matrix1, matrix2, productMatrix);
        // printMatrix(productMatrix);
    }

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}

void fillMatrices(int matrix1[][MATRIX_SIZE], int matrix2[][MATRIX_SIZE])
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            matrix1[i][j] = rand() % 10;
            matrix2[i][j] = rand() % 10;
        }
    }
}

void printMatrix(int matrix[][MATRIX_SIZE])
{
    printf("[");
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        printf("%s", i == 0 ? "[" : " [");
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("%d%s", matrix[i][j], j != MATRIX_SIZE - 1 ? ", " : "");
        }
        printf("%s", i != MATRIX_SIZE - 1 ? "]\n" : "]");
    }
    printf("]\n");
}

void multiplyMatrices(int matrix1[][MATRIX_SIZE], int matrix2[][MATRIX_SIZE], int productMatrix[][MATRIX_SIZE])
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            productMatrix[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                productMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}
