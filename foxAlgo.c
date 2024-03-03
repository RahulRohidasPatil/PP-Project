#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


#define MATRIX_SIZE 4
#define NUM_PROCESSES 16

void multiplyMatFunc(double**, double** , double**, int);

void printMatrix(double** arr,int n_bar){

    for(int i=0; i<n_bar; i++){

        for(int j=0; j<n_bar; j++){

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

void Set_to_n(double** arr, double n,int n_bar){

    for(int i=0; i<n_bar; i++){

        for(int j = 0; j<n_bar; j++){

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

int Local_matrix_allocate(double ***array, int n_bar){
    /* allocate the n*m contiguous items */
    double *p = (double *)malloc(n_bar*n_bar*sizeof(double));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (double **)malloc(n_bar*sizeof(double*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n_bar; i++)
       (*array)[i] = &(p[i*n_bar]);

    return 0;




    // double** temp = (double**)malloc(nbar*nbar*sizeof(double*));

    // for(int i=0; i<MATRIX_SIZE; i++){

    //     temp[i] = (double*)malloc(nbar*sizeof(double));   

    // }
    // return temp;
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
    Set_to_n(local_C,0,n_bar);


    //calculate address for circular shift of B
    int source = (grid->my_row + 1) % grid->q;
    int dest = (grid-> my_row + grid->q - 1) % grid->q;

    //set aside storage for the brodcast block of A
    Local_matrix_allocate(&tempA,n_bar);

    

    for(stage=0;stage < grid->q; stage++){
        bcast_root = (grid->my_row + stage) % grid->q;

        if(bcast_root == grid->my_col){
            MPI_Bcast(local_A[0],n_bar * n_bar,MPI_DOUBLE,bcast_root,grid->row_comm);
            multiplyMatFunc(local_A,local_B,local_C,n_bar);
        }else{
            MPI_Bcast(tempA[0],n_bar * n_bar,MPI_DOUBLE,bcast_root,grid->row_comm);
            multiplyMatFunc(tempA,local_B,local_C,n_bar);
        }
        
        MPI_Sendrecv_replace(local_B[0],n_bar * n_bar,MPI_DOUBLE,dest,0,source,0,grid->col_comm, &status);
    }
}

void multiplyMatFunc(double** matrix1, double** matrix2, double** productMatrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                productMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

int free2dchar(double ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

int main(int argc, char* argv[]){
    srand(time(0));

    double** matrix_a;

    double** matrix_b;

    double** matrix_c;

    int rank;

    int size;

    int ierr;

    int procgridSize = MATRIX_SIZE / (MATRIX_SIZE / sqrt(NUM_PROCESSES));
    int procElemSize = MATRIX_SIZE / sqrt(NUM_PROCESSES);



    MPI_Status status;

    MPI_Init(&argc, &argv);

    GRID_INFO_T grid;
    grid.p=NUM_PROCESSES;  
    Setup_grid(&grid);
    

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0){
        Local_matrix_allocate(&matrix_a,MATRIX_SIZE);
        Local_matrix_allocate(&matrix_b,MATRIX_SIZE);
        Local_matrix_allocate(&matrix_c,MATRIX_SIZE);

        int k = 1;
        for(int i=0; i<MATRIX_SIZE; i++){
            int temp = k;
            if(i > 0 && i % 2 == 0){
                k += 2;
                temp = k;
            }
            for(int j = 0; j<MATRIX_SIZE; j++){
                if(j > 0 && j % 2 == 0){
                    k += 1;
                }
                matrix_a[i][j] = k;

            }
            k = temp;
        }

        Set_to_n(matrix_b,1,MATRIX_SIZE);        

        // fillArray(matrix_a);
        // fillArray(matrix_b);
        // Set_to_n(matrix_a,1,MATRIX_SIZE);
        // Set_to_n(matrix_b,1,MATRIX_SIZE);

        printMatrix(matrix_a,MATRIX_SIZE);
        printMatrix(matrix_b,MATRIX_SIZE);

        // printf("matrix init completed by rank 0");
        // free(matrix_a);
        // free(matrix_b);
        // free(matrix_c);

    }

    double ** localA;
    double ** localB;
    double ** localC;
    Local_matrix_allocate(&localA,procElemSize);
    Local_matrix_allocate(&localB,procElemSize);
    Local_matrix_allocate(&localC,procElemSize);

    int sizes[2] = {MATRIX_SIZE,MATRIX_SIZE}; // Main matrix size
    int subsizes[2] = {procElemSize, procElemSize};
    int startsWith[2] = {0,0};
    MPI_Datatype type, subArrType;
    MPI_Type_create_subarray(2,sizes,subsizes,startsWith,MPI_ORDER_C,MPI_DOUBLE, &type);
    MPI_Type_create_resized(type,0,MATRIX_SIZE/procgridSize * sizeof(double),&subArrType);
    MPI_Type_commit(&subArrType);   

    double *globalptrA=NULL;
    if (rank == 0) globalptrA = &(matrix_a[0][0]);
    double *globalptrB=NULL;
    if (rank == 0) globalptrB = &(matrix_b[0][0]);
    double *globalptrC=NULL;
    if (rank == 0) globalptrC = &(matrix_c[0][0]);

    int sendcounts[procgridSize*procgridSize];
    int displs[procgridSize*procgridSize]; 
    
    if(grid.my_rank== 0){
        for (int i=0;i<procgridSize* procgridSize;i++) sendcounts[i] = 1;
        int disp = 0;
        for(int i =0; i < procgridSize;i++){
            for(int j=0;j<procgridSize;j++){
                displs[i*procgridSize+j] = disp;
                disp += 1; 
            }
            disp += ((MATRIX_SIZE/procgridSize) - 1) * procgridSize;   
        }
    }


    MPI_Scatterv(globalptrA, sendcounts,displs,subArrType,&(localA[0][0]),MATRIX_SIZE * MATRIX_SIZE / (procgridSize * procgridSize), MPI_DOUBLE,0,grid.comm);
    MPI_Scatterv(globalptrB, sendcounts,displs,subArrType,&(localB[0][0]),MATRIX_SIZE * MATRIX_SIZE / (procgridSize * procgridSize), MPI_DOUBLE,0,grid.comm);
    
    
    Fox(MATRIX_SIZE,&grid,localA,localB,localC);
    //Set_to_n(localC,5,procgridSize);

    MPI_Gatherv(&(localC[0][0]), MATRIX_SIZE* MATRIX_SIZE/(procgridSize*procgridSize),MPI_DOUBLE,globalptrC,sendcounts,displs,subArrType,0,grid.comm);

    if(rank==0){
        printMatrix(matrix_c,MATRIX_SIZE);
        Set_to_n(matrix_c,0,MATRIX_SIZE);
        multiplyMatFunc(matrix_a,matrix_b,matrix_c,MATRIX_SIZE);
        printMatrix(matrix_c,MATRIX_SIZE);
    }

    MPI_Type_free(&subArrType);

    free2dchar(&localA);
    free2dchar(&localB);
    free2dchar(&localC);

    if(rank == 0){
        free2dchar(&matrix_a);
        free2dchar(&matrix_b);
        free2dchar(&matrix_c);
    }
    

    ierr = MPI_Finalize();
}

