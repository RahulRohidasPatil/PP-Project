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



int main(int argc, char* argv[]){

}
