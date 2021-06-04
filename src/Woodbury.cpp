// Woodbury.cpp
// Woodbury 2x2 and 3x3 kernels
//
// Woodbury matrix identity:
// (S + U * V)^{-1} = S^{-1} - S^{-1} * U * B^{-1} * V * S^{-1}
// B := 1 + V * C,  2 x 2
// C := S^{-1} * U, dim x 2
#include "Woodbury.hpp"
#include "Helpers.hpp"

// Woodbury 2x2 kernel: 
void WB2(double *Slater_inv, unsigned int Dim,
         double *Updates, unsigned int *Updates_index) {
  std::cerr << "Called Woodbury 2x2 kernel" << std::endl;

  // Construct V from Updates_index
  unsigned int *V = new unsigned int[2 * Dim]{0};
  for (unsigned int i = 0; i < Dim; i++) {
    if (i == Updates_index[0] - 1) V[i] = 1;
    if (Dim + i == Updates_index[1] - 1) V[Dim + i] = 1;
  }

  // Compute B from U, V and Slater_inv
  double B[4]{0};
  double Binv[4]{0};
  double C[4]{0};
  double D[2*Dim]{0};
    // Compute C
    for (unsigned int i = 0; i < Dim; i++) {
      for (unsigned int j = 0; j < 2; j++) {
        for (unsigned int k = 0; k < Dim; k++) {
          C[i * Dim + j] += Slater_inv[i * Dim + k] * Updates[k * Dim + j];
        }
      }
    }
    // Compute B
    for (unsigned int i = 0; i < 2; i++) {
      for (unsigned int j = 0; j < 2; j++) {
        for (unsigned int k = 0; k < Dim; k++) {
          B[i * Dim + j] += (i == j) + V[i * Dim + k] * C[k * Dim + j];
        }
      }
    }

  // Invert B with explicit formula for 2x2 inversion
  double idet = 1.0 / (B[0] * B[3] - B[1] * B[2]);
  Binv[0] = idet * B[3];
  Binv[1] = -1.0 * idet * B[1];
  Binv[2] = -1.0 * idet * B[2];
  Binv[3] = idet * B[0];

  // Compute (S + U * V)^{-1} with Woobury identity
  
}

// Woodbury 3x3 kernel
void WB3(double *Slater_inv, unsigned int Dim,
         double *Updates, unsigned int *Updates_index) {
           std::cerr << "Called Woodbury 3x3 kernel" << std::endl;

           // Construct V from Updates_index

           // Compute B from U, V and Slater_inv

           // Invert B with explicit formula for 3x3 inversion

           // Compute (S + U * V)^{-1} with Woobury identity
           
         }
