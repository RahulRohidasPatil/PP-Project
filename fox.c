#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MATRIX_SIZE 4

typedef struct
{
    MPI_Comm comm;     // Communicator for entire grid
    MPI_Comm row_comm; // Communicator for my row
    MPI_Comm col_comm; // Communicator for my col
    int my_row;        // My row number
    int my_col;        // My column number
    int my_rank;       // My rank in the grid comm
} GRID_INFO_T;

void allocateMatrix(int ***matrix);
void fillMatrix(int **matrix);
void printMatrix(int **matrix);
void freeMatrix(int **matrix);
void Setup_grid(GRID_INFO_T *);
void Fox(GRID_INFO_T grid, int local_A, int local_B, int *local_C);

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int size, rank, **matrixA, **matrixB, **matrixC, localA, localB, localC = 0;
    GRID_INFO_T grid;
    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (sqrt(size) != MATRIX_SIZE)
    {
        if (rank == 0)
            printf("Number of processors must match the square of MATRIX_SIZE.\n");
        MPI_Finalize();
        return 0;
    }

    allocateMatrix(&matrixA);
    allocateMatrix(&matrixB);
    allocateMatrix(&matrixC);

    if (rank == 0)
    {
        fillMatrix(matrixA);
        printf("Matrix A:\n");
        printMatrix(matrixA);

        fillMatrix(matrixB);
        printf("Matrix B:\n");
        printMatrix(matrixB);
    }

    Setup_grid(&grid);
    MPI_Scatter(matrixA[0], 1, MPI_INT, &localA, 1, MPI_INT, 0, grid.comm);
    MPI_Scatter(matrixB[0], 1, MPI_INT, &localB, 1, MPI_INT, 0, grid.comm);

    Fox(grid, localA, localB, &localC);

    MPI_Gather(&localC, 1, MPI_INT, matrixC[0], 1, MPI_INT, 0, grid.comm);

    if (rank == 0)
    {
        printf("Matrix C:\n");
        printMatrix(matrixC);
        freeMatrix(matrixA);
        freeMatrix(matrixB);
        freeMatrix(matrixC);
    }

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}

void allocateMatrix(int ***matrix)
{
    *matrix = (int **)malloc(MATRIX_SIZE * sizeof(int *));
    int *block = (int *)malloc(MATRIX_SIZE * MATRIX_SIZE * sizeof(int));
    for (int i = 0; i < MATRIX_SIZE; i++)
        (*matrix)[i] = block + i * MATRIX_SIZE;
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
    free(matrix[0]);
    free(matrix);
}

void Setup_grid(GRID_INFO_T *grid)
{
    int dimensions[2] = {MATRIX_SIZE, MATRIX_SIZE};
    int wrapAround[2] = {1, 1};
    int coordinates[2];
    int freeCoords[2];

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrapAround, 1, &(grid->comm));
    MPI_Comm_rank(grid->comm, &(grid->my_rank));
    MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);

    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

    // Set up row communicators
    freeCoords[0] = 0;
    freeCoords[1] = 1;
    MPI_Cart_sub(grid->comm, freeCoords, &(grid->row_comm));

    // set up column communicators
    freeCoords[0] = 1;
    freeCoords[1] = 0;
    MPI_Cart_sub(grid->comm, freeCoords, &(grid->col_comm));
}

void Fox(GRID_INFO_T grid, int local_A, int local_B, int *local_C)
{
    int tempA = local_A;
    int source = (grid.my_row + 1) % MATRIX_SIZE;
    int dest = (grid.my_row + MATRIX_SIZE - 1) % MATRIX_SIZE;
    MPI_Status status;

    for (int stage = 0; stage < MATRIX_SIZE; stage++)
    {
        int bcast_root = (grid.my_row + stage) % MATRIX_SIZE;
        if (bcast_root == grid.my_col)
        {
            MPI_Bcast(&local_A, 1, MPI_INT, bcast_root, grid.row_comm);
            *local_C += local_A * local_B;
        }
        else
        {
            MPI_Bcast(&tempA, 1, MPI_INT, bcast_root, grid.row_comm);
            *local_C += tempA * local_B;
        }
        MPI_Sendrecv_replace(&local_B, 1, MPI_INT, dest, 0, source, 0, grid.col_comm, &status);
    }
}