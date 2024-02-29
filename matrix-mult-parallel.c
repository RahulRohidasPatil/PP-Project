#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_SIZE 576

void fillMatrices(int matrix1[][MATRIX_SIZE], int matrix2[][MATRIX_SIZE]);
void printMatrix(int array[][MATRIX_SIZE]);

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int size, rank, n = MATRIX_SIZE * MATRIX_SIZE;
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
    }
    MPI_Bcast(matrix2, n, MPI_INT, 0, MPI_COMM_WORLD);

    int local_n = n / size;
    int received_row_count = local_n / MATRIX_SIZE;
    int received_buffer[received_row_count][MATRIX_SIZE];
    MPI_Scatter(matrix1, local_n, MPI_INT, received_buffer, local_n, MPI_INT, 0, MPI_COMM_WORLD);

    int result_buffer[received_row_count][MATRIX_SIZE];
    for (int i = 0; i < received_row_count; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            result_buffer[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; k++)
                result_buffer[i][j] += received_buffer[i][k] * matrix2[k][j];
        }
    }

    MPI_Gather(result_buffer, local_n, MPI_INT, productMatrix, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    // if (rank == 0)
    // {
    //     printf("Product Matrix:\n");
    //     printMatrix(productMatrix);
    // }

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
