#include<stdio.h>
#include<omp.h>
long long m = 100000000;

int main(int argc, char *argv[]) {
  long long i;
  double sum;

  sum = 0.0;

  #pragma omp parallel for num_threads(10) reduction (+ : sum)
  for (i = 1; i <= m; i++) {
    sum = sum + i;
  }

  printf("sum = %.0f\n", sum);
  return 0;
}