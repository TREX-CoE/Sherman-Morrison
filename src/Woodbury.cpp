// Woodbury.cpp
// Woodbury 2x2 and 3x3 kernels
//
// Woodbury matrix identity:
// (S + U * V)^{-1} = S^{-1} - S^{-1} * U * B^{-1} * V * S^{-1}
// B := 1 + V * C,  2 x 2
// C := S^{-1} * U, dim x 2
// D := V * S^{-1}, 2 x dim
//
// All matrices are stored in row-major order

#include "Woodbury.hpp"
#include "Helpers.hpp"

// Woodbury 2x2 kernel:
bool WB2(double *Slater_inv, unsigned int Dim, double *Updates,
         unsigned int *Updates_index) {
  std::cerr << "Called Woodbury 2x2 kernel" << std::endl;

  // Construct V from Updates_index
  unsigned int V[2 * Dim]{0}; // 2 x Dim matrix stored in row-major order
  V[Updates_index[0] - 1] = 1;
  V[Dim + Updates_index[1] - 1] = 1;

  // Compute C
  double C[2 * Dim]{0};
  matMul2(Slater_inv, Updates, C, Dim, Dim, 2);
  // Compute B
  double B[4]{0};
  matMul2(V, C, B, 2, Dim, 2);
  // Compute 1 + B
  B[0] += 1;
  B[3] += 1;

  // Invert 1 + B with explicit formula for 2x2 inversion
  double idet = 1.0 / (B[0] * B[3] - B[1] * B[2]);
  double Binv[4]{0};
  Binv[0] = idet * B[3];
  Binv[1] = -1.0 * idet * B[1];
  Binv[2] = -1.0 * idet * B[2];
  Binv[3] = idet * B[0];

  // Check if determinant of inverted matrix is not zero
  double det = B[0] * B[3] - B[1] * B[2];
  if (std::fabs(det) < threshold()) {
    std::cerr << "Determinant approached 0!" << std::endl;
    return false;
  }

  // Compute (S + U * V)^{-1} with Woobury identity
  double D[2 * Dim]{0};
  matMul2(V, Slater_inv, D, 2, Dim, Dim);
  double tmp[2 * Dim]{0};
  matMul2(Binv, D, tmp, 2, 2, Dim);
  double tmp2[Dim * Dim]{0};
  matMul2(C, tmp, tmp2, Dim, 2, Dim);
  for (unsigned int i = 0; i < Dim * Dim; i++) {
    Slater_inv[i] -= tmp2[i];
  }
  return true;
}

// Woodbury 3x3 kernel
bool WB3(double *Slater_inv, unsigned int Dim, double *Updates,
         unsigned int *Updates_index) {
  std::cerr << "Called Woodbury 3x3 kernel" << std::endl;
  
  // Construct V from Updates_index
  unsigned int V[3 * Dim]{0}; // 2 x Dim matrix stored in row-major order
  V[Updates_index[0] - 1] = 1;
  V[Dim + Updates_index[1] - 1] = 1;
  V[2 * Dim + Updates_index[2] - 1] = 1;

  // Compute C
  double C[3 * Dim]{0};
  matMul2(Slater_inv, Updates, C, Dim, Dim, 3);
  // Compute B
  double B[9]{0};
  matMul2(V, C, B, 3, Dim, 3);
  // Compute 1 + B
  B[0] += 1;
  B[4] += 1;
  B[8] += 1;

  double Binv[9];
  Binv[0] = B[4] * B[8] - B[5] * B[7];
  Binv[3] = B[5] * B[6] - B[3] * B[8];
  Binv[6] = B[3] * B[7] - B[4] * B[6];

  Binv[1] = B[2] * B[7] - B[1] * B[8];
  Binv[4] = B[0] * B[8] - B[2] * B[6];
  Binv[7] = B[1] * B[6] - B[0] * B[7];

  Binv[2] = B[1] * B[5] - B[2] * B[4];
  Binv[5] = B[2] * B[3] - B[0] * B[5];
  Binv[8] = B[0] * B[4] - B[1] * B[3];

  // Check if determinant of inverted matrix is not zero
  // If so, exigt and return false.
  double det;
  det = B[0] * (B[4] * B[8] - B[5] * B[7]) -
        B[1] * (B[3] * B[8] - B[5] * B[6]) + B[2] * (B[3] * B[7] - B[4] * B[6]);
  if (std::fabs(det) < threshold()) {
    std::cerr << "Determinant approached 0!" << std::endl;
    return false;
  }

  // Compute (S + U * V)^{-1} with Woobury identity
  double D[3 * Dim]{0};
  matMul2(V, Slater_inv, D, 3, Dim, Dim);
  double tmp[3 * Dim]{0};
  matMul2(Binv, D, tmp, 3, 3, Dim);
  double tmp2[Dim * Dim]{0};
  matMul2(C, tmp, tmp2, Dim, 3, Dim);
  for (unsigned int i = 0; i < Dim * Dim; i++) {
    Slater_inv[i] -= tmp2[i];
  }
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
