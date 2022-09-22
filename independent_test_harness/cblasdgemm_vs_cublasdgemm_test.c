/*
Compile with:

nvc \
-I$NV_CUDA_MATH_PATH/11.7/include \
-L$NV_CUDA_MATH_PATH/11.7/lib64 \
-L${MKLROOT}/lib/intel64 \
-lmkl_intel_lp64 \
-lmkl_sequential \
-lmkl_core \
-lpthread \
-lm \
-ldl \
-lcublas \
-mp \
-target=gpu \
cblasdgemm_vs_cublasdgemm_test.c \
-o cblasdgemm_vs_cublasdgemm_test

*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mkl_lapacke.h>
#include <mkl.h>

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "debug.h"

int main() {

  cublasHandle_t handle;
  if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
    fprintf(stdout, "cuBLAS initialization failed!\n");
    exit(EXIT_FAILURE);
  }

  uint16_t M = 3;
  uint16_t N = 2;
  uint16_t K = 4;

  double *a = malloc(M * K * sizeof(double));  // M x K = 3 x 4 
  double *acm = malloc(M * K * sizeof(double)); // col-major stored
  double *b = malloc(K * N * sizeof(double));  // K x N = 4 x 2
  double *bcm = malloc(K * N * sizeof(double)); // col-major stored
  double *c = malloc(M * N * sizeof(double));  // M x N = 3 x 2

  a[0] = 1, a[1] = 2, a[2] = 3, a[3] = 4, a[4] = 5, a[5] = 6, a[6] = 7, a[7] = 8; a[8] = 9, a[9] = 10, a[10] = 11, a[11] = 12;
  acm[0] = 1, acm[1] = 5, acm[2] = 9, acm[3] = 2, acm[4] = 6, acm[5] = 10, acm[6] = 3, acm[7] = 7; acm[8] = 11, acm[9] = 4, acm[10] = 8, acm[11] = 12;
  b[0] = 13, b[1] = 14, b[2] = 15, b[3] = 16, b[4] = 17, b[5] = 18, b[6] = 19, b[7] = 20;
  bcm[0] = 13, bcm[1] = 15, bcm[2] = 17, bcm[3] = 19, bcm[4] = 14, bcm[5] = 16, bcm[6] = 18, bcm[7] = 20;

  uint16_t lda   = K;
  uint16_t ldacm = M;
  uint16_t ldb   = N;
  uint16_t ldbcm = K;
  uint16_t ldc   = N;

  double alpha = 1.0, beta = 0.0;
  

  cblas_dgemm(CblasRowMajor,
              CblasNoTrans, CblasNoTrans,
              M, N, K,
              alpha, a, lda, b, ldb,
              beta, c, ldc);
  print_m(c, M, N, ldc, "c_cblas_dgemm");
  

  memset(c, 0, M*N*sizeof(double));
  #pragma omp target enter data map(to:a[0:M*K], b[0:K*N], c[0:M*N])
  #pragma omp target data use_device_ptr(a, b, c)
  {
    int cublasError = cublasDgemm(handle,
                                  CUBLAS_OP_N, CUBLAS_OP_N,
                                  N, M, K,
                                  &alpha, b, ldb, a, lda,
                                  &beta, c, ldc);
  }
  #pragma omp target exit data map(from:c[0:M*N])
  print_m(c, M, N, ldc, "c_cublasDgemm");
  

  memset(c, 0, M*N*sizeof(double));
  ldc = M; // ldc : N -> M, because cublasDgemm stores result in col-maj
  #pragma omp target enter data map(to:acm[0:M*K], bcm[0:K*N], c[0:M*N])
  #pragma omp target data use_device_ptr(acm, bcm, c)
  {
    int cublasError = cublasDgemm(handle,
                                  CUBLAS_OP_N, CUBLAS_OP_N,
                                  M, N, K,
                                  &alpha, acm, ldacm, bcm, ldbcm,
                                  &beta, c, ldc);
  }
  #pragma omp target exit data map(from:c[0:M*N])
  print_m_t(c, M, N, ldc, "c_col-maj_cublasDgemm");


  free(a);
  free(acm);
  free(b);
  free(bcm);
  free(c);
  return 0;
}
