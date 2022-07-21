#include "meuk.h"
#include <stdint.h>
#include <assert.h>

void print_matrix(double *A, const uint64_t LDS, const uint64_t Dim) {
  for (uint64_t i = 0; i < LDS * Dim; i++) {
    printf("%f\n", A[i]);
  }
  printf("\n");
}

double frobenius_norm2(double *A, const uint64_t LDS, const uint64_t Dim) {
  double sum2 = 0;
  for (uint64_t i = 0; i < LDS * Dim; i++) sum2 += A[i] * A[i];
  return sum2;
}

double frobenius_norm(double *A, const uint64_t LDS, const uint64_t Dim) {
  double sum2 = frobenius_norm2(A, LDS, Dim);
  return sqrt(sum2);
}

double max_norm(double *A, const uint64_t LDS, const uint64_t Dim) {
  double largest = 0;
  for (uint64_t i = 0; i < LDS * Dim; i++) {
    double elm = A[i];
    double felm = fabs(elm);
    if (elm != elm) return -1.0; // Return a negative norm when NaN found
    if (felm > largest) largest = felm;
  }
  return largest;
}

double condition_number(double *A, double *Ainv, const uint64_t LDS, const uint64_t Dim) {
  double norm_A = frobenius_norm(A, LDS, Dim);
  double norm_Ainv = frobenius_norm(Ainv, LDS, Dim);
  return fabs(norm_A) * fabs(norm_Ainv);
}

void read_uint(hid_t file_id, const char *key, uint64_t *data) {
  herr_t rc;
  hid_t dataset_id = H5Dopen2(file_id, key, H5P_DEFAULT);
  assert(dataset_id >= 0 && "H5Dopen2");
  rc = H5Dread(dataset_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  assert(rc >= 0 && "H5Dread");
  rc = H5Dclose(dataset_id);
  assert(rc >= 0 && "H5Dclose");
}

void read_double(hid_t file_id, const char *key, double *data) {
  herr_t rc;
  hid_t  dataset_id = H5Dopen2(file_id, key, H5P_DEFAULT);
  assert(dataset_id >= 0 && "H5Dopen2");
  rc = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  assert(rc >= 0 && "H5Dread");
  rc = H5Dclose(dataset_id);
  assert(rc >= 0 && "H5Dclose");
}

void update_slater_matrix(const uint64_t LDS, const uint64_t Dim,
                          const uint64_t N_updates, const double *Updates,
                          const uint64_t *Updates_index, double *Slater) {

  for (uint32_t i = 0; i < N_updates; i++) {
    uint32_t col = Updates_index[i] - 1;
    for (uint32_t j = 0; j < Dim; j++) {
      Slater[col * Dim + j]  += Updates[i * LDS + j];
    }
  }
}

uint32_t check_error(const uint64_t LDS, const uint64_t Dim, double *Slater_invT,
                        double *Slater, const double tolerance) {

  double res[Dim*Dim];

  for (uint32_t i = 0; i < Dim; i++) {
    for (uint32_t j = 0; j < Dim; j++) {
      res[i * Dim + j] = 0;
      for (uint32_t k = 0; k < Dim; k++) {
        res[i * Dim + j] += Slater[i * Dim + k] * Slater_invT[k * LDS + j];
      }
    }
  }

  for (uint32_t i = 0; i < Dim; i++) {
    for (uint32_t j = 0; j < Dim; j++) {
      double elm = res[i * Dim + j];
      if (elm != elm) return 1; // found a NaN!
      if (i == j && fabs(elm - 1.0) > tolerance) return 1;
      if (i != j && fabs(elm) > tolerance) return 1;
    }
  }

  return 0;
}

void matmul(double *a, double *b, double *prod, const uint64_t LDS, const uint64_t Dim) {
  for (uint32_t i = 0; i < Dim; i++) {
    for (uint32_t j = 0; j < Dim; j++) {
      prod[i * Dim + j] = 0;
      for (uint32_t k = 0; k < Dim; k++) {
        prod[i * Dim + j] += a[i * Dim + k] * b[k * LDS + j];
      }
    }
  }
}

int32_t check_error_better(const double max, const double tolerance) {
  if (max < 0) return -1; // When max was a NaN
  else if (max < tolerance) return 0; // Good
  else return 1; // Too big
}

void residual(double *a, double *res, const uint64_t Dim) {
  for (uint32_t i = 0; i < Dim; i++) {
    for (uint32_t j = 0; j < Dim; j++) {
      if (i == j) res[i * Dim + j] = a[i * Dim + j] - 1.0;
      else res[i * Dim + j] = a[i * Dim + j];
    }
  }
}

uint32_t test_kernel(char *version, const uint64_t LDS, const uint64_t Dim,
                     const uint64_t N_updates, const double *Updates,
                     const uint64_t *Updates_index, const double breakdown, const double tolerance,
                     double *Slater, double *Slater_inv, double *determinant) {
  uint32_t rc = 0;
  // if (version[0] == 'a') { // Anthony
  //   const double *Upds;
  //   const uint64_t *Ui;
  //   for (int i = 0; i < LDS * Dim; i++) Slater_inv[i] *= *determinant;
  //   for (int j = 0; j < N_updates; j++) {
  //     Upds = &Updates[j * LDS];
  //     Ui = &Updates_index[j];
  //     detupd(Dim, LDS, Upds, Ui, Slater_inv, determinant);
  //     if (determinant == 0) printf("TEST_KERNEL: det_update21 failed\n");
  //   }
  //   for (int i = 0; i < LDS * Dim; i++) Slater_inv[i] /= *determinant;
  //   update_slater_matrix(LDS, Dim, N_updates, Updates, Updates_index, Slater);
  //   rc = check_error(LDS, Dim, Slater_inv, Slater, tolerance);
  //   if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  // } else if (version[0] == 'n') { // Naive
  if (version[0] == 'n') { // Naive
    rc = qmckl_sherman_morrison(LDS, Dim, N_updates, Updates, Updates_index,
                                breakdown, Slater_inv, determinant);
    if (rc != 0) printf("TEST_KERNEL: qmckl_sherman_morrison failed\n");
    update_slater_matrix(LDS, Dim, N_updates, Updates, Updates_index, Slater);
    rc = check_error(LDS, Dim, Slater_inv, Slater, tolerance);
    if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  } else if (version[0] == 's') { // Splitting
    rc = qmckl_sherman_morrison_splitting(LDS, Dim, N_updates, Updates,
                                          Updates_index, breakdown, Slater_inv,
                                          determinant);
    if (rc != 0) printf("TEST_KERNEL: qmckl_sherman_morrison_splitting failed\n");
    update_slater_matrix(LDS, Dim, N_updates, Updates, Updates_index, Slater);
    rc = check_error(LDS, Dim, Slater, Slater_inv, tolerance);
    if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  } else if (version[0] == 'b') { // Blocked
    rc = qmckl_sherman_morrison_smw32s(LDS, Dim, N_updates, Updates,
                                       Updates_index, breakdown, Slater_inv,
                                       determinant);
    if (rc != 0) printf("TEST_KERNEL: qmckl_sherman_morrison_smw32s failed\n");
    update_slater_matrix(LDS, Dim, N_updates, Updates, Updates_index, Slater);
    rc = check_error(LDS, Dim, Slater, Slater_inv, tolerance);
    if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  }
  return rc;
}
