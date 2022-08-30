// This file provide a simple API to generate random cycles
// The goal of this api is to provide a way for other programs to generate random cycles at runtime
// instead of requiring a pre-generated dataset

// Don't use "pragma once" to ensure compatibility with any compiler
#ifndef C_API_HEADER_GUARD
#define C_API_HEADER_GUARD

// This file is included from both C++ and C, so ensure everything here has C linkage
// if need be
#ifdef __cplusplus
#include <cstdint>
extern "C" {
#else
#include <stdint.h>
#endif


typedef struct {
  // Total number of updates in the cycle
  uint32_t n_update;
  // Dim is the number of rows
  uint32_t dim;
  // lds is the number of columns
  uint32_t lds;
  double *slater_inverse_t;
  double *slater_matrix;
  // Matrix containing all the updates !The updates are ADDITIVE!
  double *updates;
  // Determinant of the slater matrix
  double determinant;
  // Indices of the updates
  uint32_t *col_update_index;
} Cycle;

Cycle *generateRandomCycle(uint32_t dim, uint32_t lds, uint32_t n_update, uint32_t n_splits,
                           double determinant_threshold, double splitting_update_noise_magnitude);

void freeCycle(Cycle **cycle);


// To match the previous opening bracket
#ifdef __cplusplus
}
#endif

#endif
