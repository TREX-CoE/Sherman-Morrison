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
bool WB2(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index, const double breakdown) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/
#ifdef DEBUG1
  std::cerr << "Called Woodbury 2x2 kernel" << std::endl;
#endif

  const unsigned int row1 = (Updates_index[0] - 1);
  const unsigned int row2 = (Updates_index[1] - 1);

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE
  // OF LAYOUT OF 'Updates' !!
  double C[2 * Dim];
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < 2; j++) {
      C[i * 2 + j] = 0;
      for (unsigned int k = 0; k < Dim; k++) {
        C[i * 2 + j] += Slater_inv[i * Dim + k] * Updates[Dim * j + k];
      }
    }
  }

  // Compute B = 1 + V * C
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (std::fabs(det) < breakdown) {
#ifdef DEBUG1
    std::cerr << "Determinant too close to zero! No inverse found."
              << std::endl;
    std::cerr << "Determinant = " << det << std::endl;
#endif
    return false;
  }

  // Compute B^{-1} with explicit formula for 2x2 inversion
  double Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // Compute tmp = B^{-1} x (V.S^{-1})
  double tmp[2 * Dim];
  for (unsigned int i = 0; i < 2; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      tmp[i * Dim + j] = Binv[i * 2] * Slater_inv[row1 * Dim + j];
      tmp[i * Dim + j] += Binv[i * 2 + 1] * Slater_inv[row2 * Dim + j];
    }
  }

  // Compute (S + U V)^{-1} = S^{-1} - C x tmp
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      Slater_inv[i * Dim + j] -= C[i * 2] * tmp[j];
      Slater_inv[i * Dim + j] -= C[i * 2 + 1] * tmp[Dim + j];
    }
  }

  return true;
}

// Woodbury 3x3 kernel
bool WB3(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index, const double breakdown) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/
#ifdef DEBUG1
  std::cerr << "Called Woodbury 3x3 kernel" << std::endl;
#endif

  const unsigned int row1 = (Updates_index[0] - 1);
  const unsigned int row2 = (Updates_index[1] - 1);
  const unsigned int row3 = (Updates_index[2] - 1);

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE
  // OF LAYOUT OF 'Updates' !!
  double C[3 * Dim];
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      C[i * 3 + j] = 0;
      for (unsigned int k = 0; k < Dim; k++) {
        C[i * 3 + j] += Slater_inv[i * Dim + k] * Updates[Dim * j + k];
      }
    }
  }

#ifdef DEBUG2
  showMatrix2(C, Dim, 3, "C = S_inv * U");
  showMatrix2(D, 3, Dim, "D = V * S_inv");
#endif

  // Compute B = 1 + V.C
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

#ifdef DEBUG2
  showMatrix2(B, 3, 3, "B = 1 + V * C");
#endif

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
#ifdef DEBUG2
  std::cerr << "Determinant of B = " << det << std::endl;
#endif
  if (std::fabs(det) < breakdown) {
#ifdef DEBUG1
    std::cerr << "Determinant too close to zero! No inverse found."
              << std::endl;
    std::cerr << "Determinant = " << det << std::endl;
#endif
    return false;
  }

  // Compute B^{-1} with explicit formula for 3x3 inversion
  double Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

#ifdef DEBUG2
  std::cerr << "Conditioning number of B = " << condition1(B, Binv, 3)
            << std::endl;
  showMatrix2(Binv, 3, 3, "Binv");
#endif

  // Compute tmp = B^{-1} x (V.S^{-1})
  double tmp[3 * Dim];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      tmp[i * Dim + j] = Binv[i * 3] * Slater_inv[row1 * Dim + j];
      tmp[i * Dim + j] += Binv[i * 3 + 1] * Slater_inv[row2 * Dim + j];
      tmp[i * Dim + j] += Binv[i * 3 + 2] * Slater_inv[row3 * Dim + j];
    }
  }

#ifdef DEBUG2
  showMatrix2(tmp, 3, Dim, "tmp = Binv * D");
#endif

  // Compute (S + U V)^{-1} = S^{-1} - C x tmp
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      Slater_inv[i * Dim + j] -= C[i * 3] * tmp[j];
      Slater_inv[i * Dim + j] -= C[i * 3 + 1] * tmp[Dim + j];
      Slater_inv[i * Dim + j] -= C[i * 3 + 2] * tmp[2 * Dim + j];
    }
  }

#ifdef DEBUG2
  showMatrix2(Slater_inv, Dim, Dim, "Slater_inv AFTER update");
#endif
  return true;
}

extern "C" {
bool WB2_f(double **linSlater_inv, unsigned int *Dim, double **linUpdates,
           unsigned int **Updates_index, const double breakdown) {
  bool ok;
  ok = WB2(*linSlater_inv, *Dim, *linUpdates, *Updates_index, breakdown);
  return ok;
}
bool WB3_f(double **linSlater_inv, unsigned int *Dim, double **linUpdates,
           unsigned int **Updates_index, const double breakdown) {
  bool ok;
  ok = WB3(*linSlater_inv, *Dim, *linUpdates, *Updates_index, breakdown);
  return ok;
}
}
