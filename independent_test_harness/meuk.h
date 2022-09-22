#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5/serial/hdf5.h>
#include "kernels.h"

typedef struct Error {
  uint32_t rc;
  uint64_t error;
} Error;

void matmul(double *a, double *b, double *prod, const uint64_t Lds, const uint64_t Dim);
void residual(double *a, double *res, const uint64_t Dim);
double frobenius_norm2(double *A, const uint64_t Lds, const uint64_t Dim);
double frobenius_norm(double *A, const uint64_t Lds, const uint64_t Dim);
double max_norm(double *A, const uint64_t Lds, const uint64_t Dim);
double condition_number(double *A, double *Ainv, const uint64_t Lds, const uint64_t Dim);
void read_uint(hid_t file_id, const char *key, uint64_t *data);
void read_double(hid_t file_id, const char *key, double *data);

static __inline__ uint64_t rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}

void update_slater_matrix(const uint64_t Lds, const uint64_t Dim,
                          const uint64_t N_updates, const double *Updates,
                          const uint64_t *Updates_index, double *Slater);

uint32_t  check_error(const uint64_t Lds, const uint64_t Dim, double *Slater_invT,
                        double *Slater, const double tolerance);

int32_t check_error_better(const double max, const double tolerance);

uint32_t test_kernel(char *version, const uint64_t Lds, const uint64_t Dim,
                     const uint64_t N_updates, const double *Updates,
                     const uint64_t *Updates_index, const double breakdown, const double tolerance,
                     double *Slater, double *Slater_inv, double *determinant);
