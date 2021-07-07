/*  Woodbury.cpp
    Woodbury 2x2 and 3x3 kernels

    Woodbury matrix identity:
    (S + U V)^{-1}  =  S^{-1} - C B^{-1} D

    All matrices are stored in row-major order
*/

#include "Woodbury.hpp"
#include "Helpers.hpp"

// #define DEBUG1
// #define DEBUG2

// Woodbury 2x2 kernel
bool WB2(double *Slater_inv, unsigned int Dim, double *Updates,
         unsigned int *Updates_index) {
  /*
      C := S^{-1} * U,    dim x 2
      B := 1 + V * C,     2 x 2
      D := V * S^{-1},    2 x dim
  */
#ifdef DEBUG1
  std::cerr << "Called Woodbury 2x2 kernel" << std::endl;
#endif

  // Compute D = V.S^{-1}
  double D[2 * Dim];
  unsigned int row1, row2;
  row1 = Updates_index[0] - 1;
  row2 = Updates_index[1] - 1;
  for (unsigned int i = 0; i < Dim; i++) {
    D[i] = Slater_inv[row1 * Dim + i];
    D[Dim + i] = Slater_inv[row2 * Dim + i];
  }

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE 
  // OF LAYOUT OF 'Updates' !!
  double C[2 * Dim];
  for(unsigned int i = 0; i < Dim; i++) {
    for(unsigned int j = 0; j < 2; j++) {
      C[i * 2 + j] = 0;
      for(unsigned int k = 0; k < Dim; k++) {
        C[i * 2 + j] += Slater_inv[i * Dim + k] * Updates[Dim * j + k];
      }
    }
  }

  // Compute B = 1 + V * C
  double B[4];
  B[0] = C[row1 * 2];
  B[1] = C[row1 * 2 + 1];
  B[2] = C[row2 * 2];
  B[3] = C[row2 * 2 + 1];
  B[0] += 1, B[3] += 1;

  // Check if determinant of inverted matrix is not zero
  double det = B[0] * B[3] - B[1] * B[2];
  if (std::fabs(det) < threshold()) {
    std::cerr << "Determinant too close to zero! No inverse found." << std::endl;
#ifdef DEBUG1
    std::cerr << "Determinant = " << det << std::endl;
#endif
    return false;
  }

  // Compute B^{-1} with explicit formula for 2x2 inversion
  double Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B[3];
  Binv[1] = -1.0 * idet * B[1];
  Binv[2] = -1.0 * idet * B[2];
  Binv[3] = idet * B[0];

  // Compute B^{-1} x D
  double tmp[2 * Dim];
  matMul2(Binv, D, tmp, 2, 2, Dim);

  // Compute C x B^{-1} x D
  double tmp2[Dim * Dim];
  matMul2(C, tmp, tmp2, Dim, 2, Dim);

  // Compute (S + U V)^{-1} = S^{-1} - C B^{-1} D
  for (unsigned int i = 0; i < Dim * Dim; i++) {
    Slater_inv[i] -= tmp2[i];
  }

  return true;
}

// Woodbury 3x3 kernel
bool  WB3(double *Slater_inv, unsigned int Dim, double *Updates,
         unsigned int *Updates_index) {
  /*
      C := S^{-1} * U,    dim x 3
      B := 1 + V * C,     3 x 3
      D := V * S^{-1},    3 x dim
  */
#ifdef DEBUG1
  std::cerr << "Called Woodbury 3x3 kernel" << std::endl;
#endif
#ifdef DEBUG2
  showMatrix2(Slater_inv, Dim, Dim, "Slater_inv BEFORE update");
  showMatrix2(Updates, 3, Dim, "Updates");
  showMatrix2(Updates_index, 1, 3, "Updates_index");
#endif

  // Compute D = V * S^{-1}
  double D[3 * Dim];
  unsigned int row1, row2, row3;
  row1 = Updates_index[0] - 1;
  row2 = Updates_index[1] - 1;
  row3 = Updates_index[2] - 1;
  for (unsigned int i = 0; i < Dim; i++) {
    D[i] = Slater_inv[row1 * Dim + i];
    D[Dim + i] = Slater_inv[row2 * Dim + i];
    D[2 * Dim + i] = Slater_inv[row3 * Dim + i];
  }

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE 
  // OF LAYOUT OF 'Updates' !!
  double C[3 * Dim];
  for(unsigned int i = 0; i < Dim; i++) {
    for(unsigned int j = 0; j < 3; j++) {
      C[i * 3 + j] = 0;
      for(unsigned int k = 0; k < Dim; k++) {
        C[i * 3 + j] += Slater_inv[i * Dim + k] * Updates[Dim * j + k];
      }
    }
  }

#ifdef DEBUG2
  showMatrix2(C, Dim, 3, "C = S_inv * U");
#endif

#ifdef DEBUG2
  showMatrix2(D, 3, Dim, "D = V * S_inv");
#endif

  // Compute B = 1 + V.C
  double B[9];
  B[0] = C[row1 * 3];
  B[1] = C[row1 * 3 + 1];
  B[2] = C[row1 * 3 + 2];
  B[3] = C[row2 * 3];
  B[4] = C[row2 * 3 + 1];
  B[5] = C[row2 * 3 + 2];
  B[6] = C[row3 * 3];
  B[7] = C[row3 * 3 + 1];
  B[8] = C[row3 * 3 + 2];
  B[0] += 1, B[4] += 1, B[8] += 1;

#ifdef DEBUG2
  showMatrix2(B, 3, 3, "B = 1 + V * C");
#endif

  // Check if determinant of B is not too close to zero
  double det;
  det = B[0] * (B[4] * B[8] - B[5] * B[7]) -
        B[1] * (B[3] * B[8] - B[5] * B[6]) +
        B[2] * (B[3] * B[7] - B[4] * B[6]);
#ifdef DEBUG2 
  std::cerr << "Determinant of B = " << det << std::endl;
#endif
  if (std::fabs(det) < threshold()) {
  // if (std::fabs(det) < 1000000) {
    std::cerr << "Determinant too close to zero! No inverse found." << std::endl;
#ifdef DEBUG1
    std::cerr << "Determinant = " << det << std::endl;
#endif
    return false;
  }

  // Compute B^{-1} with explicit formula for 3x3 inversion
  double Binv[9], idet = 1.0 / det;
  Binv[0] =   ( B[4] * B[8] - B[7] * B[5] ) * idet;
  Binv[1] = - ( B[1] * B[8] - B[7] * B[2] ) * idet;
  Binv[2] =   ( B[1] * B[5] - B[4] * B[2] ) * idet;
  Binv[3] = - ( B[3] * B[8] - B[6] * B[5] ) * idet;
  Binv[4] =   ( B[0] * B[8] - B[6] * B[2] ) * idet;
  Binv[5] = - ( B[0] * B[5] - B[3] * B[2] ) * idet;
  Binv[6] =   ( B[3] * B[7] - B[6] * B[4] ) * idet;
  Binv[7] = - ( B[0] * B[7] - B[6] * B[1] ) * idet;
  Binv[8] =   ( B[0] * B[4] - B[3] * B[1] ) * idet;

#ifdef DEBUG2
  std::cerr << "Conditioning number of B = " << condition1(B, Binv, 3) << std::endl;
  showMatrix2(Binv, 3, 3, "Binv");
#endif

  // Compute B^{-1} x D
  double tmp[3 * Dim];
  matMul2(Binv, D, tmp, 3, 3, Dim);

#ifdef DEBUG2
  showMatrix2(tmp, 3, Dim, "tmp = Binv * D");
#endif

  // Compute C x B^{-1} x D
  double tmp2[Dim * Dim];
  matMul2(C, tmp, tmp2, Dim, 3, Dim);

#ifdef DEBUG2
  showMatrix2(tmp2, Dim, Dim, "tmp2 = C * tmp");
#endif

  // Compute (S + U V)^{-1} = S^{-1} - C B^{-1} D
  for (unsigned int i = 0; i < Dim * Dim; i++) {
    Slater_inv[i] -= tmp2[i];
  }

#ifdef DEBUG2
  showMatrix2(Slater_inv, Dim, Dim, "Slater_inv AFTER update");
#endif
  return true;
}

extern "C" {
bool WB2_f(double **linSlater_inv, unsigned int *Dim, double **linUpdates,
           unsigned int **Updates_index) {
  WB2(*linSlater_inv, *Dim, *linUpdates, *Updates_index);
}
bool WB3_f(double **linSlater_inv, unsigned int *Dim, double **linUpdates,
           unsigned int **Updates_index) {
  WB3(*linSlater_inv, *Dim, *linUpdates, *Updates_index);
}
}
