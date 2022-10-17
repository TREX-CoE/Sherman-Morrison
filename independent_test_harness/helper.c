#include "helper.h"
#include <stdint.h>
#include <assert.h>

#ifdef HAVE_CUBLAS_OFFLOAD
  cublasHandle_t init_cublas() {
    cublasHandle_t handle;
    if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
      fprintf(stdout, "cuBLAS initialization failed!\n");
      exit(EXIT_FAILURE);
    }
    return handle;
  }

  cusolverDnHandle_t init_cusolver() {
    cusolverDnHandle_t handle;
    if (cusolverDnCreate(&handle) != CUSOLVER_STATUS_SUCCESS) {
      fprintf(stdout, "cuSOLVER initialization failed!\n");
      exit(EXIT_FAILURE);
    }
    return handle;
  }
#endif

void copy(double* Slater_invT_copy, uint64_t Lds, double* tmp, uint64_t Dim) {
  for (uint32_t i = 0; i < Dim; i++) {
    for (uint32_t j = 0; j < Lds; j++) {
      if (j < Dim) Slater_invT_copy[i * Lds + j] = tmp[i * Dim + j];
      else Slater_invT_copy[i * Lds + j] = 0.0;
    }
  }
}

void update(double* slaterT,double* upds, uint64_t* ui, uint64_t nupds,uint64_t Dim, u_int64_t Lds) {
  for (int i = 0; i < nupds; i++) {
    int col = ui[i] - 1;
    for (int j = 0; j < Dim; j++) {
      slaterT[col + j * Dim] += upds[i * Lds + j];
    }
  }
}

void convert(double* upds, uint64_t nupds, uint64_t* ui, double* slaterT, uint64_t Dim, u_int64_t Lds) {
  for (int i = 0; i < nupds; i++) {
    int col = ui[i] - 1;
    for (int j = 0; j < Lds; j++) {
      upds[i * Lds + j] -= slaterT[col + j * Dim];
    }
  }
}

void transpose(double* a, uint16_t lda, double *b, uint16_t ldb, uint16_t m, uint16_t n)
{
  for(uint16_t i = 0; i < m; i++)
  {
    for( uint16_t j = 0; j < n; j++)
    {
      b[j * ldb + i] = a[i * lda + j];
    }
  }
}

double get_determinant(uint32_t cycle, hid_t file_id) {
  char det_key[32];
  sprintf(det_key, "/cycle_%d/determinant", cycle);
  double determinant;
  read_double(file_id, det_key, &determinant);
  return determinant;
}

double* get_slater_inv(uint32_t cycle, hid_t file_id, uint64_t Dim, u_int64_t Lds) {
  char slater_inv_key[32];
  sprintf(slater_inv_key, "/cycle_%d/slater_inverse_t", cycle);
  double *slater_inv = malloc(sizeof *slater_inv * Dim * Lds);
  read_double(file_id, slater_inv_key, slater_inv);
  return slater_inv;
}

double* get_slater(uint32_t cycle, hid_t file_id, uint64_t Dim, u_int64_t Lds) {
  char slater_key[32];
  sprintf(slater_key, "/cycle_%d/slater_matrix", cycle);
  double *slater = malloc(sizeof *slater * Dim * Lds);
  read_double(file_id, slater_key, slater);
  return slater;
}

double* get_upds(uint32_t cycle, hid_t file_id, uint64_t nupds, u_int64_t Lds) {
  char upds_key[32];
  sprintf(upds_key, "/cycle_%d/updates", cycle);
  double *upds = malloc(sizeof *upds * Lds * nupds);
  read_double(file_id, upds_key, upds);
  return upds;
}

uint64_t* get_upd_idcs(uint32_t cycle, hid_t file_id, uint64_t nupds) {
  char upd_idx_key[32];
  sprintf(upd_idx_key, "/cycle_%d/col_update_index", cycle);
  uint64_t* uis = malloc(sizeof *uis * nupds);
  read_uint(file_id, upd_idx_key, uis);
  return uis;
}

uint64_t get_dim(uint32_t cycle, hid_t file_id) {
  char dim_key[32];
  sprintf(dim_key, "/cycle_%d/slater_matrix_dim", cycle);
  uint64_t Dim;
  read_uint(file_id, dim_key, &Dim);
  return Dim;
}

uint64_t get_nupdates(uint32_t cycle, hid_t file_id) {
  char nupds_key[32];
  sprintf(nupds_key, "/cycle_%d/nupdates", cycle);
  uint64_t N_updates;
  read_uint(file_id, nupds_key, &N_updates);
  return N_updates;
}

double frobenius_norm2(double *A, const uint64_t Lds, const uint64_t Dim) {
  double sum2 = 0;
  for (uint64_t i = 0; i < Lds * Dim; i++) sum2 += A[i] * A[i];
  return sum2;
}

double frobenius_norm(double *A, const uint64_t Lds, const uint64_t Dim) {
  double sum2 = frobenius_norm2(A, Lds, Dim);
  return sqrt(sum2);
}

double max_norm(double *A, const uint64_t Lds, const uint64_t Dim) {
  double largest = 0;
  for (uint64_t i = 0; i < Lds * Dim; i++) {
    double elm = A[i];
    double felm = fabs(elm);
    if (elm != elm) return -1.0; // Return a negative norm when NaN found
    if (felm > largest) largest = felm;
  }
  return largest;
}

double condition_number(double *A, double *Ainv, const uint64_t Lds, const uint64_t Dim) {
  double norm_A = frobenius_norm(A, Lds, Dim);
  double norm_Ainv = frobenius_norm(Ainv, Lds, Dim);
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

void update_slater_matrix(const uint64_t Lds, const uint64_t Dim,
                          const uint64_t N_updates, const double *Updates,
                          const uint64_t *Updates_index, double *Slater) {

  for (uint32_t i = 0; i < N_updates; i++) {
    uint32_t col = Updates_index[i] - 1;
    for (uint32_t j = 0; j < Dim; j++) {
      Slater[col * Dim + j]  += Updates[i * Lds + j];
    }
  }
}

uint32_t check_error(const uint64_t Lds, const uint64_t Dim, double *Slater_invT,
                        double *Slater, const double tolerance) {

  double* res = malloc(sizeof *res * Dim * Dim);
  double alpha = 1.0, beta = 0.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              Dim, Dim, Dim,
              alpha, Slater, Dim, Slater_invT, Lds,
              beta, res, Dim);

  for (uint32_t i = 0; i < Dim; i++) {
    for (uint32_t j = 0; j < Dim; j++) {
      double elm = res[i * Dim + j];
      if (elm != elm) return 1; // found a NaN!
      if (i == j && fabs(elm - 1.0) > tolerance) return 1;
      if (i != j && fabs(elm) > tolerance) return 1;
    }
  }
  free(res);
  return 0;
}

//void matmul(double *a, double *b, double *prod, const uint64_t Lds, const uint64_t Dim) {
//  for (uint32_t i = 0; i < Dim; i++) {
//    for (uint32_t j = 0; j < Dim; j++) {
//      prod[i * Dim + j] = 0;
//      for (uint32_t k = 0; k < Dim; k++) {
//        prod[i * Dim + j] += a[i * Dim + k] * b[k * Lds + j];
//      }
//    }
//  }
//}

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

uint32_t test_kernel(char *version, const uint64_t Lds, const uint64_t Dim,
                     const uint64_t N_updates, const double *Updates,
                     const uint64_t *Updates_index, const double breakdown, const double tolerance,
                     double *Slater, double *Slater_inv, double *determinant) {
  uint32_t rc = 0;

  if (version[0] == 'n') { // Naive
    rc = qmckl_sherman_morrison(Lds, Dim, N_updates, Updates, Updates_index,
                                breakdown, Slater_inv, determinant);
    if (rc != 0) printf("TEST_KERNEL: qmckl_sherman_morrison failed\n");
    update_slater_matrix(Lds, Dim, N_updates, Updates, Updates_index, Slater);
    rc = check_error(Lds, Dim, Slater_inv, Slater, tolerance);
    if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  } else if (version[0] == 's') { // Splitting
    rc = qmckl_sherman_morrison_splitting(Lds, Dim, N_updates, Updates,
                                          Updates_index, breakdown, Slater_inv,
                                          determinant);
    if (rc != 0) printf("TEST_KERNEL: qmckl_sherman_morrison_splitting failed\n");
    update_slater_matrix(Lds, Dim, N_updates, Updates, Updates_index, Slater);
    rc = check_error(Lds, Dim, Slater, Slater_inv, tolerance);
    if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  } else if (version[0] == 'b') { // Blocked
    rc = qmckl_sherman_morrison_smw32s(Lds, Dim, N_updates, Updates,
                                       Updates_index, breakdown, Slater_inv,
                                       determinant);
    if (rc != 0) printf("TEST_KERNEL: qmckl_sherman_morrison_smw32s failed\n");
    update_slater_matrix(Lds, Dim, N_updates, Updates, Updates_index, Slater);
    rc = check_error(Lds, Dim, Slater, Slater_inv, tolerance);
    if (rc != 0) printf("TEST_KERNEL: check_error failed\n");
  }
  return rc;
}
