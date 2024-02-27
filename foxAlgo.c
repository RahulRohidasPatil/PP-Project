#include<mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define MATRIX_SIZE 8000
#define NUM_PROCESSES 10

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
    MPI_Cart_sub(grid->comm,free_coords, &(grid->rom_comm));

    //set up col communicators
    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid->comm,free_coords, &(grid->col_comm));
}


void Fox(
    int         n,
    GRID_INFO_T* grid,
    LOCAL_MATRIX_T* local_A,
    LOCAL_MATRIX_T* local_B,
    LOCAL_MATRIX_T* local_C
){

    LOCAL_MATIX_T* tempA;
    int stage;
    int bcast_root;
    int n_bar;
    int score;
    int dist;
    MPI_Status status;


    n_bar = n / grid->q;
    Set_to_zero(local_C);


    //calculate address for circular shift of B
    source = (grid->my_row + 1) % grid->q;
    dest = (grid-> my_row + grid->q - 1) % grid->q;

    //set aside storage for the brodcast block of A
    temp_A = Local_matrix_allocate(n_bar);

    for(stage=0;stage < grid->q; stage++){
        bcast_root = (grid->my_row + stage) % grid->q;

        if(bcast_root == grid->my_col){
            MPI_Bcast(local_A,1,local_matrix_mpi_t,bcast_root,grid->row_comm);
            Local_matrix_multiply(local_A,local_B,local_C);
        }else{
            MPI_Bcast(temp_A,1,local_matrix_mpi_t,bcast_root,grid->row_comm);
            Local_matrix_multiply(local_A,local_B,local_C);
        }

        MPI_Sendrecv_replace(local_B,1,local_matrix_mpi_t,dest,0,source,0,grid->col_comm, &status);

    }
}

int main(int argc, char* argv[]){

}
