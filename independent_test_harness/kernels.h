#include <mkl_lapacke.h>
#include <mkl.h>

#define HAVE_CUBLAS_OFFLOAD

#ifdef HAVE_CUBLAS_OFFLOAD
#include <stdio.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#endif

lapack_int inverse(double *A, uint64_t Dim, uint64_t LDS);

int min(int a, int b);

uint32_t qmckl_sherman_morrison(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant);

uint32_t qmckl_sherman_morrison_splitting(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant);

uint32_t qmckl_sherman_morrison_smw32s(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant);

uint32_t qmckl_woodbury_3(const uint64_t vLDS, const uint64_t vDim,
                          const double *__restrict __attribute__((aligned(8)))
                          Updates,
                          const uint64_t *__restrict Updates_index,
                          const double breakdown,
                          double *__restrict __attribute__((aligned(8)))
                          Slater_inv,
                          double *__restrict determinant);

uint32_t qmckl_woodbury_k(const uint64_t vLDS,
                          const uint64_t vDim,
                          const uint64_t N_updates,
                          const double *__restrict __attribute__((aligned(8))) Updates,
                          const uint64_t *__restrict Updates_index,
                          const double breakdown,
                          double *__restrict __attribute__((aligned(8))) Slater_inv,
                          double *__restrict determinant);

#ifdef HAVE_CUBLAS_OFFLOAD
uint32_t qmckl_woodbury_k_cublas_offload(const uint64_t vLDS,
          const uint64_t vDim,
          const uint64_t N_updates,
          const double *__restrict __attribute__((aligned(8))) Updates,
          const uint64_t *__restrict Updates_index,
          const double breakdown,
          double *__restrict __attribute__((aligned(8))) Slater_inv,
          double *__restrict determinant);
#endif

uint32_t qmckl_woodbury_2(const uint64_t vLDS, const uint64_t vDim,
                          const double *__restrict __attribute__((aligned(8)))
                          Updates,
                          const uint64_t *__restrict Updates_index,
                          const double breakdown,
                          double *__restrict __attribute__((aligned(8)))
                          Slater_inv,
                          double *__restrict determinant);

void detupd(const uint64_t Dim, const uint64_t LDS,
            const double *__restrict __attribute__((aligned(8))) Updates,
            const uint64_t *__restrict Updates_index, 
            double *__restrict __attribute__((aligned(8))) Slater_inv,
            double *__restrict determinant);

uint32_t qmckl_sherman_morrison_later(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant);
            