#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


#define MATRIX_SIZE 4
#define NUM_PROCESSES 2

void multiplyMatFunc(double**, double** , double**);

void printMatrix(double** arr){

    for(int i=0; i<MATRIX_SIZE; i++){

        for(int j=0; j<MATRIX_SIZE; j++){

            printf("%.2f ", arr[i][j]);

        }

        printf("\n");

    }

    printf("\n");

}

void fillArray(double** arr){

    for(int i=0; i<MATRIX_SIZE; i++){

        for(int j = 0; j<MATRIX_SIZE; j++){

            int random = rand();

            arr[i][j] = (double)random/RAND_MAX;

        }

    }

}

void Set_to_n(double** arr, double n){

    for(int i=0; i<MATRIX_SIZE; i++){

        for(int j = 0; j<MATRIX_SIZE; j++){

            arr[i][j] = n;

        }

    }
}

typedef struct {
    int     p; // Total number of processes
    MPI_Comm comm; // Communicator for entire grid
    MPI_Comm row_comm; // Communicator for my row
    MPI_Comm col_comm; // Communicator for my col
    int     q; //Order of grid
    int     my_row; // My row number
    int     my_col; // My col number
    int     my_rank; // My rank in the grid comm
} GRID_INFO_T;

void Setup_grid(
    GRID_INFO_T* grid /* out */ 
){

    int old_rank;
    int dimensions[2];
    int wrap_arounds[2];
    int coordinates[2];
    int free_coords[2];

    //Set up global grid information

    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
    MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

    // we assume p is a perfect square
    grid->q = (int) sqrt((double) grid->p);
    dimensions[0] = dimensions[1] = grid->q;

    // we want a circular shift in second dimension 
    // Dont care about first 
    wrap_arounds[0] = wrap_arounds[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_arounds,1,&(grid->comm));
    MPI_Comm_rank(grid->comm,&(grid->my_rank));
    MPI_Cart_coords(grid->comm,grid->my_rank,2,coordinates);
    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

    //set up row communicators
    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(grid->comm,free_coords, &(grid->row_comm));

    //set up col communicators
    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid->comm,free_coords, &(grid->col_comm));
}

double** Local_matrix_allocate(int nbar){
    double** temp = (double**)malloc(nbar*nbar*sizeof(double*));

    for(int i=0; i<MATRIX_SIZE; i++){

        temp[i] = (double*)malloc(MATRIX_SIZE*sizeof(double));   

    }
    return temp;
}


void Fox(
    int         n,
    GRID_INFO_T* grid,
    double** local_A,
    double** local_B,
    double** local_C
){

    double** tempA;
    int stage;
    int bcast_root;
    int n_bar;
    int score;
    int dist;
    MPI_Status status;


    n_bar = n / grid->q;
    Set_to_n(local_C,0);


    //calculate address for circular shift of B
    int source = (grid->my_row + 1) % grid->q;
    int dest = (grid-> my_row + grid->q - 1) % grid->q;

    //set aside storage for the brodcast block of A
    tempA = Local_matrix_allocate(n_bar);

    for(stage=0;stage < grid->q; stage++){
        bcast_root = (grid->my_row + stage) % grid->q;

        if(bcast_root == grid->my_col){
            MPI_Bcast(local_A,1,MPI_DOUBLE,bcast_root,grid->row_comm);
            multiplyMatFunc(local_A,local_B,local_C);
        }else{
            MPI_Bcast(tempA,1,MPI_DOUBLE,bcast_root,grid->row_comm);
            multiplyMatFunc(local_A,local_B,local_C);
        }

        MPI_Sendrecv_replace(local_B,1,MPI_DOUBLE,dest,0,source,0,grid->col_comm, &status);

    }
}

void multiplyMatFunc(double** matrix1, double** matrix2, double** productMatrix)
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

int main(int argc, char* argv[]){
    srand(time(0));

    double** matrix_a;

    double** matrix_b;

    double** matrix_c;

    int rank;

    int size;

    int ierr;

    MPI_Status status;

    MPI_Init(&argc, &argv);

    GRID_INFO_T grid;
    grid.p=NUM_PROCESSES;    
    
    Setup_grid(&grid);
    

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("my rank %d : row = %d, col = %d\n",grid.my_rank,grid.my_row,grid.my_col);
    

    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0){
        matrix_a = Local_matrix_allocate(MATRIX_SIZE);
        matrix_b = Local_matrix_allocate(MATRIX_SIZE);
        matrix_c = Local_matrix_allocate(MATRIX_SIZE);

        fillArray(matrix_a);
        fillArray(matrix_b);

        printMatrix(matrix_a);
        printMatrix(matrix_b);

        printf("matrix init completed by rank 0");
        free(matrix_a);
        free(matrix_b);
        free(matrix_c);

    }else{
        printf("rank %d ready to serve",rank);
    }

    ierr = MPI_Finalize();
}

