#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int size, rank;
  /* No MPI calls before this */
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("Hello World from process %d of %d!\n", rank, size);

  MPI_Finalize();
  /* No MPI calls after this */
  return 0;
}