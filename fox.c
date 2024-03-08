#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct
{
    MPI_Comm comm;     // Communicator for entire grid
    MPI_Comm row_comm; // Communicator for my row
    MPI_Comm col_comm; // Communicator for my col
    int q;             // Order of grid
    int my_row;        // My row number
    int my_col;        // My column number
    int my_rank;       // My rank in the grid comm
} GRID_INFO_T;

int grid_order, local_matrix_size, local_matrix_elements;

void setup_grid(GRID_INFO_T *);
void allocate_matrix(int ***, int);
void fill_matrix(int **, int);
void print_matrix(int **, int);
void fox(GRID_INFO_T, int **, int **, int **);
void set_to_zero(int **);
void local_matrix_multiply(int **, int **, int **);
void freeMatrix(int **);

int main(int argc, char *argv[])
{
    // Declaration of required variables
    GRID_INFO_T grid;
    int size, rank, matrix_size;
    int **matrix_A, **matrix_B, **matrix_C, **local_A, **local_B, **local_C;

    /* No MPI calls before this */
    MPI_Init(&argc, &argv);

    // Get number of processes and current process rank within MPI communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Throw error if number of processes is not a perfect square
    grid_order = sqrt(size);
    if (grid_order * grid_order != size)
    {
        if (rank == 0)
            printf("Number of processes must be a perfect square\n");
        MPI_Finalize();
        return 0;
    }

    // Seed the random number generator with current system time
    srand(time(NULL));

    // Initialize MPI grid communicator
    setup_grid(&grid);

    // Set size of matrix using command line argument
    matrix_size = atoi(argv[1]);

    // Allocate Memory to Matrices
    allocate_matrix(&matrix_A, matrix_size);
    allocate_matrix(&matrix_B, matrix_size);
    allocate_matrix(&matrix_C, matrix_size);

    // Fill Random data in both matrices
    if (grid.my_rank == 0)
    {
        fill_matrix(matrix_A, matrix_size);
        printf("Matrix A:\n");
        print_matrix(matrix_A, matrix_size);

        fill_matrix(matrix_B, matrix_size);
        printf("Matrix B:\n");
        print_matrix(matrix_B, matrix_size);
    }

    // Allocate memory to submatrices
    local_matrix_size = matrix_size / grid.q;
    local_matrix_elements = local_matrix_size * local_matrix_size;
    allocate_matrix(&local_A, local_matrix_size);
    allocate_matrix(&local_B, local_matrix_size);
    allocate_matrix(&local_C, local_matrix_size);

    // Create a custom MPI Type to scatter submatrices
    int sizes[2] = {matrix_size, matrix_size};
    int subsizes[2] = {local_matrix_size, local_matrix_size};
    int starts[2] = {0, 0};
    MPI_Datatype local_matrix_mpi_t;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &local_matrix_mpi_t);
    MPI_Type_create_resized(local_matrix_mpi_t, 0, local_matrix_size * sizeof(int), &local_matrix_mpi_t);
    MPI_Type_commit(&local_matrix_mpi_t);

    // Scatter submatrices to all processes
    int sendcounts[size];
    int displs[size];
    if (grid.my_rank == 0)
    {
        int displ = local_matrix_size * grid.q;
        for (int i = 0; i < grid.q; i++)
            for (int j = 0; j < grid.q; j++)
            {
                sendcounts[grid.q * i + j] = 1;
                displs[grid.q * i + j] = displ * i + j;
            }
    }
    MPI_Scatterv(*matrix_A, sendcounts, displs, local_matrix_mpi_t, *local_A, local_matrix_elements, MPI_INT, 0, grid.comm);
    MPI_Scatterv(*matrix_B, sendcounts, displs, local_matrix_mpi_t, *local_B, local_matrix_elements, MPI_INT, 0, grid.comm);

    // Start Fox Algorithm
    fox(grid, local_A, local_B, local_C);

    // Gather all the results
    MPI_Gatherv(*local_C, local_matrix_elements, MPI_INT, *matrix_C, sendcounts, displs, local_matrix_mpi_t, 0, grid.comm);

    // Print the result matrix
    if (grid.my_rank == 0)
    {
        printf("Matrix C:\n");
        print_matrix(matrix_C, matrix_size);
    }

    // Free all the allocated memory
    freeMatrix(matrix_A);
    freeMatrix(matrix_B);
    freeMatrix(matrix_C);
    freeMatrix(local_A);
    freeMatrix(local_B);
    freeMatrix(local_C);
    MPI_Type_free(&local_matrix_mpi_t);

    MPI_Finalize();
    /* No MPI calls after this */
    return 0;
}

void setup_grid(GRID_INFO_T *grid)
{
    // Declaration of required variables
    int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    int free_coords[2];

    // Initialize required variables
    grid->q = grid_order;
    dimensions[0] = dimensions[1] = grid->q;

    // We want a circular shift in second dimension
    // Don't care about first
    wrap_around[0] = wrap_around[1] = 1;

    // Create a 2D Cartesian MPI Communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &grid->comm);
    MPI_Comm_rank(grid->comm, &grid->my_rank);
    MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);

    // Extract row and column indices from process coordinates
    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

    // Set up row communicators
    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(grid->comm, free_coords, &grid->row_comm);

    // set up column communicators
    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid->comm, free_coords, &grid->col_comm);
}

void allocate_matrix(int ***matrix, int size)
{
    *matrix = malloc(size * sizeof(int *));
    int *block = malloc(size * size * sizeof(int));
    for (int i = 0; i < size; i++)
        (*matrix)[i] = block + i * size;
}

void fill_matrix(int **matrix, int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[i][j] = rand() % 10;
}

void print_matrix(int **matrix, int size)
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

void fox(GRID_INFO_T grid, int **local_A, int **local_B, int **local_C)
{
    // Storage for the submatrix of A used during the current stage
    int **tempA;

    // Fill Local Result Matrix with zero
    set_to_zero(local_C);

    int source = (grid.my_row + 1) % grid.q;
    int dest = (grid.my_row + grid.q - 1) % grid.q;

    // Set aside storage for the broadcast block of A
    allocate_matrix(&tempA, local_matrix_size);

    for (int stage = 0; stage < grid.q; stage++)
    {
        int bcast_root = (grid.my_row + stage) % grid.q;
        if (bcast_root == grid.my_col)
        {
            MPI_Bcast(*local_A, local_matrix_elements, MPI_INT, bcast_root, grid.row_comm);
            local_matrix_multiply(local_A, local_B, local_C);
        }
        else
        {
            MPI_Bcast(*tempA, local_matrix_elements, MPI_INT, bcast_root, grid.row_comm);
            local_matrix_multiply(tempA, local_B, local_C);
        }
        MPI_Sendrecv_replace(*local_B, local_matrix_elements, MPI_INT, dest, 0, source, 0, grid.col_comm, MPI_STATUS_IGNORE);
    }
}

void set_to_zero(int **matrix)
{
    for (int i = 0; i < local_matrix_size; i++)
        for (int j = 0; j < local_matrix_size; j++)
            matrix[i][j] = 0;
}

void local_matrix_multiply(int **local_A, int **local_B, int **local_C)
{
    for (int i = 0; i < local_matrix_size; i++)
        for (int j = 0; j < local_matrix_size; j++)
            for (int k = 0; k < local_matrix_size; k++)
                local_C[i][j] += local_A[i][k] * local_B[k][j];
}

void freeMatrix(int **matrix)
{
    free(*matrix);
    free(matrix);
}
