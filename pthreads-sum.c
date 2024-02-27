#include <pthread.h>
#include <stdio.h>


#define NUM_OF_THREADS 10

const long long n = 100000000;

typedef struct {
  // in
  long long start;
  long long end;
  // out
  double partSolution;
} input_t;

void *sumFn(void *_data) {
  double sum = 0.0;
  input_t *inout = (input_t *)_data;
  for (int i = inout->start; i < inout->end; i++) {
    sum = sum + (double)i;
  }
  inout->partSolution = sum;
}

int main(int argc, char *argv[]) {
  pthread_t thread[NUM_OF_THREADS];  // Array of threads
  input_t
      in[NUM_OF_THREADS];  // create for every thread an array entry with its
                           // own struct to avoid overriding by other thread
  long stepWidth = n / NUM_OF_THREADS;
  for (int j = 0; j < NUM_OF_THREADS; j++) {
    in[j].start = (j * stepWidth) + 1;
    in[j].end = ((j + 1) * stepWidth) + 1;
    pthread_create(&thread[j], NULL, sumFn, &in[j]);
  }

  for (int j = 0; j < NUM_OF_THREADS; j++) pthread_join(thread[j], NULL);

  double sum = 0.0;
  for (int j = 0; j < NUM_OF_THREADS; j++) {
    printf("partSolution for %d : %.0f\n",j,in[j].partSolution);
    sum += in[j].partSolution;
  }

  printf("sum = %.0f\n", sum);
  return 0;
}