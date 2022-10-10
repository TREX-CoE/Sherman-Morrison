#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include "kernels.h"
//#include "debug.h"

extern uint64_t n_splits;
extern uint64_t block_fail;
extern uint64_t recursive_calls;

int min(int a, int b) {
  return (a > b) ? b : a;
}

uint32_t qmckl_sherman_morrison(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  double __attribute__((aligned(8))) C[Dim];
  double __attribute__((aligned(8))) D[Lds];

  uint32_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x u_l
    for (uint32_t i = 0; i < Dim; i++) {
      C[i] = 0.0;
#pragma ivdep
#pragma vector aligned
      for (uint32_t j = 0; j < Lds; j++) {
        C[i] += Slater_inv[i * Lds + j] * Updates[l * Lds + j]; // regular mat-vec product, but actually working on S_inv^T * U_l.
      }
    }

    // Denominator: v_l^T * C
    const int cui = Updates_index[l] - 1;
    double den = 1.0 + C[cui];

    if (fabs(den) < breakdown) {
      return 1;
    }
    double iden = 1.0 / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

#pragma ivdep
#pragma vector aligned
    for (uint32_t j = 0; j < Lds; j++) {
      D[j] = Slater_inv[cui * Lds + j]; // selecting proper column of v_l^T * S_inv
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint32_t i = 0; i < Dim; i++) {
#pragma ivdep
#pragma vector aligned
      for (uint32_t j = 0; j < Lds; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * Lds + j] -= update;
      }
    }
    l += 1;
  }
  return 0;
}

/*
COMPUTE S^{-1}P - CB^{-1}D  : Dim x LDS,
where S^{-1}P               : Dim x LDS,
      C := S^{-1}PP^TU      : Dim x 2,
      B := 1 + VC           : 2 x 2,
      D := VS^{-1}P         : 2 x LDS,
      P^TU                  : LDS x 2,
      V                     : 2 x Dim
*/
uint32_t qmckl_woodbury_2(const uint64_t vLDS, const uint64_t vDim,
                          const double *__restrict __attribute__((aligned(8)))
                          Updates,
                          const uint64_t *__restrict Updates_index,
                          const double breakdown,
                          double *__restrict __attribute__((aligned(8)))
                          Slater_inv,
                          double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  const uint32_t row1 = (Updates_index[0] - 1);
  const uint32_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : Dim x 2
  double __attribute__((aligned(8))) C[2 * Dim];
  for (uint32_t i = 0; i < Dim; i++) {
    C[i * 2] = 0;
    C[i * 2 + 1] = 0;
#pragma ivdep
#pragma vector aligned
    for (uint32_t k = 0; k < Lds; k++) {
      C[i * 2]     += Slater_inv[i * Lds + k] * Updates[k];
      C[i * 2 + 1] += Slater_inv[i * Lds + k] * Updates[Lds + k];
    }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
    return 1;
  }

  // Update det(S) when passed
  if (determinant != NULL)
    *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x LDS
  double __attribute__((aligned(8))) tmp[2 * Lds];
  double *__restrict r1dim = &(Slater_inv[row1 * Lds]);
  double *__restrict r2dim = &(Slater_inv[row2 * Lds]);
#pragma ivdep
#pragma vector aligned
  for (uint32_t j = 0; j < Lds; j++) {
    tmp[j]       = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
    tmp[Lds + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : Dim x Lds
  for (uint32_t i = 0; i < Dim; i++) {
#pragma ivdep
#pragma vector aligned
    for (uint32_t j = 0; j < Lds; j++) {
      Slater_inv[i * Lds + j] -= C[i * 2]     * tmp[j];
      Slater_inv[i * Lds + j] -= C[i * 2 + 1] * tmp[Lds + j];
    }
  }

  return 0;
}

/*
COMPUTE (S^T)^{-1} - CB^{-1}D : Dim x LDS,
where S^T                     : Dim x LDS,
      C := (S^T)^{-1}U        : Dim x 3,
      B := 1 + VC             : 3 x 3,
      D := V(S^T)^{-1}        : 3 x LDS,
      U                       : LDS x 3,
      V                       : 3 x Dim
*/
uint32_t  qmckl_woodbury_3(const uint64_t vLDS, const uint64_t vDim,
                           const double *__restrict __attribute__((aligned(8)))
                           Updates,
                           const uint64_t *__restrict Updates_index,
                           const double breakdown,
                           double *__restrict __attribute__((aligned(8)))
                           Slater_inv,
                           double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  const uint32_t row1 = (Updates_index[0] - 1);
  const uint32_t row2 = (Updates_index[1] - 1);
  const uint32_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : Dim x 3
  double __attribute__((aligned(8))) C[3 * Dim];
  for (uint32_t i = 0; i < Dim; i++) {
    C[i * 3] = 0;
    C[i * 3 + 1] = 0;
    C[i * 3 + 2] = 0;
#pragma ivdep
#pragma vector aligned
    for (uint32_t k = 0; k < Lds; k++) {
      C[i * 3] += Slater_inv[i * Lds + k] * Updates[k];
      C[i * 3 + 1] += Slater_inv[i * Lds + k] * Updates[Lds + k];
      C[i * 3 + 2] += Slater_inv[i * Lds + k] * Updates[2 * Lds + k];
    }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
    return 1;
  }

  // Update det(Slater) if passed
  if (determinant != NULL)
    *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] =  (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] =  (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] =  (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] =  (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] =  (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x LDS
  double __attribute__((aligned(8))) tmp[3 * Lds];
  double *__restrict r1dim = &(Slater_inv[row1 * Lds]);
  double *__restrict r2dim = &(Slater_inv[row2 * Lds]);
  double *__restrict r3dim = &(Slater_inv[row3 * Lds]);
#pragma ivdep
#pragma vector aligned
  for (uint32_t j = 0; j < Lds; j++) {
    tmp[j]           = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
    tmp[Lds + j]     = Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
    tmp[2 * Lds + j] = Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : Dim x Lds
  for (uint32_t i = 0; i < Dim; i++) {
#pragma ivdep
#pragma vector aligned
    for (uint32_t j = 0; j < Lds; j++) {
      Slater_inv[i * Lds + j] -= C[i * 3]     * tmp[j];
      Slater_inv[i * Lds + j] -= C[i * 3 + 1] * tmp[Lds + j];
      Slater_inv[i * Lds + j] -= C[i * 3 + 2] * tmp[2 * Lds + j];
    }
  }

  return 0;
}

/*
COMPUTE S^{-1} - C B^{-1} D : Dim x LDS,
where S^{-1}                : Dim x LDS,
      C := S^{-1} U         : Dim x K, dgemm
      B := 1 + V C          : K x K, copy
      D := V S^{-1}         : K x LDS, copy
      U                     : LDS x K,
      V                     : K x Dim
      tmp := B^{-1} D       : K x LDS, dgemm
      S = S - C tmp         : Dim x LDS, dgemm
*/
uint32_t qmckl_woodbury_k(const uint64_t vLDS,
                          const uint64_t vDim,
                          const uint64_t N_updates,
                          const double *__restrict __attribute__((aligned(8))) Updates,
                          const uint64_t *__restrict Updates_index,
                          const double breakdown,
                          double *__restrict __attribute__((aligned(8))) Slater_inv,
                          double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  // Compute C = S^{-1} U : Dim x K : standard dgemm
  double *C = calloc(1, Dim * N_updates * sizeof(double));
  double alpha = 1.0, beta = 0.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
              Dim, N_updates, Lds,
              alpha, Slater_inv, Lds, Updates, Lds,
              beta, C, N_updates);

  // Construct B = 1 + V C : K x K, construct D = V S^{-1} : K x LDS
  double* B = calloc(1, sizeof *B * N_updates * N_updates);
  double* D = calloc(1, sizeof *D * N_updates * Lds);
  for (uint32_t i = 0; i < N_updates; i++) {
    const uint32_t row = Updates_index[i] - 1;
    for (uint32_t j = 0; j < N_updates  ; j++) B[i * N_updates + j] = C[row * N_updates + j] + (i == j);
    for (uint32_t j = 0; j < Lds; j++) D[i * Lds + j] = Slater_inv[row * Lds + j];
  }

  // Compute determinant by LU decomposition
  int* pivot  = calloc(1, sizeof *pivot * N_updates);
  (void) LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N_updates, N_updates, B, N_updates, pivot);

  bool swap = false; uint32_t j = 0; double det = 1.0f;
  for (uint32_t i = 0; i < N_updates; i++) {
    swap = (bool)(pivot[i] - (i + 1));     // swap = {0->false: no swap, >0->true: swap}
    j += (uint32_t)swap;                   // count # of swaps
    det *= B[i * (N_updates + 1)];         // prod. of diag elm. of B
  }
  if (fabs(det) < breakdown) return 1;  // check if determinant of B is too close to zero. If so, exit early.
  if (determinant) {                       // update det(Slater) if determinant != NULL
    if ((j & 1) != 0) det = -det;          // multiply det with -1 if # of swaps is odd
    *determinant *= det;
  }

  // Compute B^{-1} with explicit formula for K x K inversion
  (void) LAPACKE_dgetri(LAPACK_ROW_MAJOR, N_updates, B, N_updates, pivot);

  // tmp1 = B^{-1} D : KxLDS = KxK X KxLDS : standard dgemm
  double* tmp1 = calloc(1, sizeof *tmp1 * N_updates * Lds);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              N_updates, Lds, N_updates,
              alpha, B, N_updates, D, Lds,
              beta, tmp1, Lds);

  // Compute S^{-1} - C * tmp1 : Dim x LDS : standard dgemm
  alpha = -1.0, beta = 1.0;
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              Dim, Lds, N_updates,
              alpha, C, N_updates, tmp1, Lds,
              beta, Slater_inv, Lds);

  free(C);
  free(B);
  free(D);
  free(tmp1);
  free(pivot);
  return 0;
}

#ifdef HAVE_CUBLAS_OFFLOAD
uint32_t qmckl_woodbury_k_cublas_offload(cublasHandle_t b_handle, cusolverDnHandle_t s_handle,
                                         const uint64_t vLDS,
                                         const uint64_t vDim,
                                         const uint64_t N_updates,
                                         const double* Updates,
                                         const uint64_t* Updates_index,
                                         const double breakdown,
                                         double* Slater_inv,
                                         double* determinant)
{
  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  bool swap;
  uint32_t j;
  double alpha, beta, det;
  int* pivot   = malloc(sizeof *pivot * N_updates);
  double* C    = malloc(sizeof *C * Dim * N_updates);
  double* B    = malloc(sizeof *B * N_updates * N_updates);
  double* Binv = malloc(sizeof *Binv * N_updates * N_updates);
  double* D    = malloc(sizeof *D * N_updates * Lds);
  double* T1   = malloc(sizeof *T1 * N_updates * Lds);
  double* T2   = malloc(sizeof *T2 * Dim * Lds);

  int workspace_size = 0, *info = NULL;
  double* workspace = NULL;
  cusolverDnDgetrf_bufferSize(s_handle, N_updates, N_updates, B, N_updates, &workspace_size);
  printf("SIZE OF CUSOLVER WORKSPACE: %d doubles of %lu byte = %lu byte\n", workspace_size, sizeof *workspace, sizeof *workspace * workspace_size);
  workspace = malloc(sizeof *workspace * workspace_size);

  #pragma omp target enter data map(to: Updates[0:Lds*N_updates], \
                                        Updates_index[0:N_updates], \
                                        Slater_inv[0:Dim*Lds])

  // Compute C <- S^{-1} U : Dim x K : standard dgemm
  alpha = 1.0f, beta = 0.0f;
  #pragma omp target enter data map(alloc: C[0:Dim*N_updates])
  #pragma omp target data use_device_ptr(Slater_inv, Updates, C)
  {
    (void) cublasDgemm_v2(b_handle,
                          CUBLAS_OP_T, CUBLAS_OP_N,
                          N_updates, Dim, Lds,
                          &alpha, Updates, Lds, Slater_inv, Lds,
                          &beta, C, N_updates);
  }
  #pragma omp target exit data map(delete: Updates[0:Lds*N_updates])

  // Construct B <- 1 + V C : K x K
  #pragma omp target enter data map(alloc: B[0:N_updates*N_updates])
  #pragma omp target teams distribute parallel for
  for (uint32_t i = 0; i < N_updates; i++) {
    const uint32_t row = Updates_index[i] - 1;
    for (uint32_t j = 0; j < N_updates  ; j++) {
      B[j * N_updates + i] = C[row * N_updates + j] + (i == j); // B NEEDS TO BE IN COL-MAJ FOR cusolverDnDgetrf !
    }
  }

  // Compute det(B) via LU(B)
  #pragma omp target enter data map(alloc: workspace[0:workspace_size], pivot[0:N_updates])
  #pragma omp target data use_device_ptr(B, workspace, pivot)
  {
    (void) cusolverDnDgetrf(s_handle, N_updates, N_updates, B, N_updates, workspace, pivot, info);
  }
  #pragma omp target exit data map(delete: workspace[0:workspace_size])
  swap = false; j = 0; det = 1.0f;
  #pragma omp target teams distribute parallel for reduction(+: j) reduction(*: det)
  for (uint32_t i = 0; i < N_updates; i++) {
    swap = (bool)(pivot[i] - (i + 1));     // swap = {0->false: no swap, >0->true: swap}
    j += (uint32_t)swap;                   // count # of swaps
    det *= B[i * (N_updates + 1)];         // prod. of diag elm. of B
  }
  if (fabs(det) < breakdown) return 1;  // check if determinant of B is too close to zero. If so, exit early.
  if (determinant) {                       // update det(Slater) if determinant != NULL
    if ((j & 1) != 0) det = -det;          // multiply det with -1 if # of swaps is odd
    *determinant *= det;
  }
  
  // Compute B^{-1} : initialise as I for solving BX=I
  #pragma omp target enter data map(alloc: Binv[0:N_updates*N_updates])
  #pragma omp target teams distribute parallel for collapse(2)
  for (int i = 0; i < N_updates; ++i) {
    for (int j = 0; j < N_updates; ++j) {
      Binv[i * N_updates + j] = (i == j);
    }
  }
  #pragma omp target data use_device_ptr(B, pivot, Binv)
  {
    (void) cusolverDnDgetrs(s_handle, CUBLAS_OP_N, N_updates, N_updates, B, N_updates, pivot, Binv, N_updates, info);
  }
  #pragma omp target exit data map(delete: B[0:N_updates*N_updates], pivot[0:N_updates])

  // Construct D = V S^{-1} : K x LDS
  #pragma omp target enter data map(alloc: D[0:N_updates*Lds])
  #pragma omp target teams distribute parallel for
  for (uint32_t i = 0; i < N_updates; i++) {
    const uint32_t row = Updates_index[i] - 1;
    for (uint32_t j = 0; j < Lds; j++) {
      D[i * Lds + j] = Slater_inv[row * Lds + j];
    }
  }
  #pragma omp target exit data map(delete: Updates_index[0:N_updates])

  // T1 <- B^{-1} D : KxLDS : standard dgemm
  #pragma omp target enter data map(alloc: T1[0:N_updates*Lds])
  #pragma omp target data use_device_ptr(D, Binv, T1)
  {
    (void) cublasDgemm_v2(b_handle,
                          CUBLAS_OP_N,
                          CUBLAS_OP_T, // REMEMBER THIS IS Binv TRANSPOSED  because of cusolverDnDgetrs CALL ON l.434 !!!
                          Lds, N_updates, N_updates,
                          &alpha, D, Lds, Binv, N_updates,
                          &beta, T1, Lds);
  }
  #pragma omp target exit data map(delete: D[0:N_updates*Lds], Binv[0:N_updates*N_updates])

  // Compute T2 <- C * T1 : Dim x LDS : standard dgemm
  #pragma omp target enter data map(alloc: T2[0:Dim*Lds])
  #pragma omp target data use_device_ptr(T1, C, T2)
  {
    (void) cublasDgemm_v2(b_handle,
                          CUBLAS_OP_N, CUBLAS_OP_N,
                          Dim, Lds, N_updates,
                          &alpha, T1, Lds, C, N_updates,
                          &beta, T2, Lds);
  }
  #pragma omp target exit data map(delete: T1[0:N_updates*Lds], C[0:Dim*N_updates])

  // Compute S^{-1} <- S^{-1} - T2 : Dim x LDS
  #pragma omp target teams distribute parallel for
  for (uint32_t i = 0; i < Dim * Lds; i++) {
    Slater_inv[i] = Slater_inv[i] - T2[i];
  }

  #pragma omp target update from(Slater_inv[0:Dim*Lds])
  #pragma omp target exit data map(delete: Slater_inv[0:Dim*Lds], T2[0:Dim*Lds])

  free(pivot);
  free(B);
  free(Binv);
  free(C);
  free(D);
  free(T1);
  free(T2);
  return 0;
}
#endif

uint32_t qmckl_slagel_splitting(
    const uint64_t vLDS, const uint64_t vDim, uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict __attribute__((aligned(8))) later_updates,
    uint64_t *__restrict later_index, uint64_t *__restrict later,
    double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  double __attribute__((aligned(8))) C[Lds];
  double __attribute__((aligned(8))) D[Lds];

  uint32_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint32_t i = 0; i < Dim; i++) {
      C[i] = 0.0;
#pragma ivdep
#pragma vector aligned
      for (uint32_t j = 0; j < Lds; j++) {
        C[i] += Slater_inv[i * Lds + j] * Updates[l * Lds + j]; // regular mat-vec product, but actually working on S_inv^T * U_l.
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0 + C[cui];
    // printf("test breakdown = %f, den = %f, C[cui] = %f, cui = %d\n", breakdown, fabs(den), C[cui], cui);
    if (fabs(den) < breakdown) { // Here is decided to split the update, or not.
      // printf("Split! breakdown = %f\n", breakdown);
      n_splits += 1;

      // U_l = U_l / 2: split the update in 2 equal halves and save the second halve
      // in later_updates
#pragma ivdep
#pragma vector aligned
      for (uint32_t i = 0; i < Lds; i++) {
        later_updates[*later * Lds + i] = Updates[l * Lds + i] / 2.0;
        C[i] /= 2.0;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0 + C[cui];
    } // From here onwards we continue with applying the first halve of the update to Slater_inv
    double iden = 1.0 / den;

    if (determinant) *determinant *= den;

    // D = v^T x S^{-1} : 1 x LDS
#pragma ivdep
#pragma vector aligned
    for (uint32_t j = 0; j < Lds; j++) {
      D[j] = Slater_inv[cui * Lds + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint32_t i = 0; i < Dim; i++) {
#pragma ivdep
#pragma vector aligned
      for (uint32_t j = 0; j < Lds; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * Lds + j] -= update;
      }
    }
    l += 1;
  }

  return 0;
}

uint32_t qmckl_sherman_morrison_splitting(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  double __attribute__((aligned(8))) later_updates[Lds * N_updates];
  uint64_t later_index[N_updates];
  uint64_t later = 0;
  // uint32_t rc;

  (void) qmckl_slagel_splitting(Lds, Dim, N_updates, Updates, Updates_index,
                                breakdown, Slater_inv, later_updates, later_index,
                                &later, determinant);
  // rc = qmckl_slagel_splitting(Lds, Dim, N_updates, Updates, Updates_index,
  //                             breakdown, Slater_inv, later_updates, later_index,
  //                             &later, determinant);
  // if (rc != 0) printf("Something when catastrophically wrong in QMCKL_SLAGEL_SPLITTING\n");

  if (later > 0) {
    recursive_calls++;
    // printf("Later > 0\n");
    (void) qmckl_sherman_morrison_splitting(Lds, Dim, later, later_updates,
                                            later_index, breakdown, Slater_inv,
                                            determinant);

    // rc = qmckl_sherman_morrison_splitting(Lds, Dim, later, later_updates,
    //                                       later_index, breakdown, Slater_inv,
    //                                       determinant);
    // if (rc != 0) printf("Something when catastrophically wrong in QMCKL_SHERMAN_MORRISON_SPLITTING\n");
  }

  return 0;
}

uint32_t qmckl_sherman_morrison_smw32s(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  double __attribute__((aligned(8))) later_updates[Lds * N_updates];
  uint64_t later_index[N_updates];
  uint64_t later = 0;
  uint32_t rc;

  if (N_updates == 4) { // Special case for 4 rank-1 updates: 2+2
    rc = qmckl_woodbury_2(Lds, Dim, Updates, Updates_index,
                          breakdown, Slater_inv, determinant);
    if (rc != 0) { // Send the entire block to slagel_splitting
      block_fail += 1;
      uint64_t l = 0;
      rc = qmckl_slagel_splitting(Lds, Dim, 2, Updates,
                                  Updates_index, breakdown, Slater_inv,
                                  later_updates + (Lds * later),
                                  later_index + later, &l, determinant);
      later += l;
    }
    rc = qmckl_woodbury_2(Lds, Dim, &Updates[2*Lds], &Updates_index[2],
                          breakdown, Slater_inv, determinant);
    if (rc != 0) { // Send the entire block to slagel_splitting
      block_fail += 1;
      uint64_t l = 0;
      rc = qmckl_slagel_splitting(Lds, Dim, 2, &Updates[2*Lds],
                                  &Updates_index[2], breakdown, Slater_inv,
                                  later_updates + (Lds * later),
                                  later_index + later, &l, determinant);
      later += l;
    }
    if (later > 0) {
      recursive_calls++;
      rc = qmckl_sherman_morrison_splitting(Lds, Dim, later, later_updates,
                                            later_index, breakdown, Slater_inv,
                                            determinant);
    }
    return 0;
  }

  // And for the other cases != 4, 6
  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with
  // Woodbury 3x3 kernel
  uint32_t n_of_3blocks = N_updates / 3;
  uint32_t remainder = N_updates % 3;
  uint32_t length_3block = 3 * Lds;

  if (n_of_3blocks > 0) {
    for (uint32_t i = 0; i < n_of_3blocks; i++) {
      const double *Updates_3block = &Updates[i * length_3block];
      const uint64_t *Updates_index_3block = &Updates_index[i * 3];
      rc = qmckl_woodbury_3(Lds, Dim, Updates_3block, Updates_index_3block,
                            breakdown, Slater_inv, determinant);
      if (rc != 0) { // Send the entire block to slagel_splitting
        // printf("QMCKL_WOODBURY_3 failed. Sending to QMCKL_SLAGEL_SPLITTING\n");
        block_fail += 1;
        uint64_t l = 0;
        rc = qmckl_slagel_splitting(Lds, Dim, 3, Updates_3block,
                                    Updates_index_3block, breakdown, Slater_inv,
                                    later_updates + (Lds * later),
                                    later_index + later, &l, determinant);
        // if (rc != 0) printf("Something when catastrophically wrong in QMCKL_SLAGEL_SPLITTING\n");
        later += l;
      }
    }
  }

  // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
  if (remainder == 2) {
    const double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    const uint64_t *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    rc = qmckl_woodbury_2(Lds, Dim, Updates_2block, Updates_index_2block,
                          breakdown, Slater_inv, determinant);
    if (rc != 0) { // Send the entire block to slagel_splitting
      // printf("QMCKL_WOODBURY_2 failed. Sending to QMCKL_SLAGEL_SPLITTING\n");
      block_fail += 1;
      uint64_t l = 0;
      rc = qmckl_slagel_splitting(Lds, Dim, 2, Updates_2block,
                                  Updates_index_2block, breakdown, Slater_inv,
                                  later_updates + (Lds * later),
                                  later_index + later, &l, determinant);
      // if (rc != 0) printf("Something when catastrophically wrong in QMCKL_SLAGEL_SPLITTING\n");
      later += l;
    }
  }

  // Apply last remaining update with slagel_splitting
  if (remainder == 1) {
    // // printf("Sending single update to QMCKL_SLAGEL_SPLITTING\n");
    const double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    const uint64_t *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    uint64_t l = 0;
    rc = qmckl_slagel_splitting(Lds, Dim, 1, Updates_1block,
                                Updates_index_1block, breakdown, Slater_inv,
                                later_updates + (Lds * later),
                                later_index + later, &l, determinant);
    // if (rc != 0) printf("Something when catastrophically wrong in QMCKL_SLAGEL_SPLITTING\n");
    later += l;
  }

  if (later > 0) {
    recursive_calls++;
    // printf("Sending remaining updates to QMCKL_SHERMAN_MORRISON_SPLITTING\n");
    rc = qmckl_sherman_morrison_splitting(Lds, Dim, later, later_updates,
                                          later_index, breakdown, Slater_inv,
                                          determinant);
    // if (rc != 0) printf("Something when catastrophically wrong in QMCKL_SHERMAN_MORRISON_SPLITTING\n");
  }
  return 0;
}

// Sherman Morrison, leaving zero denominators for later
uint32_t qmckl_sherman_morrison_later(
    const uint64_t vLDS, const uint64_t vDim, const uint64_t N_updates,
    const double *__restrict __attribute__((aligned(8))) Updates,
    const uint64_t *__restrict Updates_index, const double breakdown,
    double *__restrict __attribute__((aligned(8))) Slater_inv,
    double *__restrict determinant) {

  const uint32_t Dim = vDim;
  const uint32_t Lds = vLDS;

  double __attribute__((aligned(8))) C[Dim];
  double __attribute__((aligned(8))) D[Lds];

  double __attribute__((aligned(8))) later_updates[Lds * N_updates];
  uint64_t later_index[N_updates];
  uint64_t later = 0;

  uint32_t l = 0;
  // For each update
  while (l < N_updates) {

    // C = A^{-1} x U_l
    for (uint32_t i = 0; i < Dim; i++) {
      C[i] = 0.0;
#pragma ivdep
#pragma vector aligned
      for (uint32_t j = 0; j < Lds; j++) {
        C[i] += Slater_inv[i * Lds + j] * Updates[l * Lds + j]; // regular mat-vec product, but actually working on S_inv^T * U_l.
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0 + C[cui];
    if (fabs(den) < breakdown) {
#pragma ivdep
#pragma vector aligned
      // for (uint32_t i = 0; i < Dim; i++) {
      for (uint32_t i = 0; i < Lds; i++) {
        later_updates[later * Lds + i] = Updates[l * Lds + i];
      }
      later_index[later] = Updates_index[l];
      later++;
      l += 1;
      continue;
    }
    double iden = 1.0 / den;

    if (determinant) *determinant *= den;

    // D = v^T x A^{-1}
#pragma ivdep
#pragma vector aligned
    for (uint32_t j = 0; j < Lds; j++) {
      D[j] = Slater_inv[cui * Lds + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint32_t i = 0; i < Dim; i++) {
#pragma ivdep
#pragma vector aligned
      for (uint32_t j = 0; j < Lds; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * Lds + j] -= update;
      }
    }
    l += 1;
  }

  if (later == N_updates) { // If all the updates have failed, exit early with an error
    return 1;
  }
  else if (later > 0) { // If some have failed, make a recursive call
    recursive_calls++;
    (void) qmckl_sherman_morrison_later(Lds, Dim, later, later_updates,
                                        later_index, breakdown, Slater_inv, determinant);
  }

  return 0;

}

// Inplace inverse n x n matrix A.
// returns:
//   ret = 0 on success
//   ret < 0 illegal argument value
//   ret > 0 singular matrix
lapack_int inverse(double *a, uint64_t m, uint64_t n) {
  int pivot[m + 1];
  lapack_int ret;
  ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, a, n, pivot);
  if (ret != 0) return ret;
  ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, pivot);
  return ret;
}
