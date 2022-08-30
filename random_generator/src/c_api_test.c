
#include "CApi.h"
#include "stdio.h"
#include "stdlib.h"

int main() {
  Cycle *cycle = generateRandomCycle(3, 5, 3, 2, 1e-3, 1e-3);
  printf("%d\n", cycle->n_update);
  printf("%d\n", cycle->dim);
  printf("%d\n", cycle->lds);
  printf("%f\n", cycle->determinant);
  // Matrix:
  for (int i = 0; i < cycle->dim; i++) {
    for (int j = 0; j < cycle->lds; j++) {
      printf("%f ", cycle->slater_matrix[i * cycle->lds + j]);
    }
    printf("\n");
  }
  printf("\n");

  // update_array:
  for (int i = 0; i < cycle->n_update; i++) {
    printf("Update %i:\n", i);
    for (int j = 0; j < cycle->lds; j++) { printf("%lf ", cycle->updates[i * cycle->lds + j]); }
    printf("\n");
  }
  freeCycle(&cycle);
  return 0;
}