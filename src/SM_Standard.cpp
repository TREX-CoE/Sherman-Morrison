// SM-Standard.cpp
// Standard Sherman Morrison with multiple updates
#include "SM_Standard.hpp"
#include "Helpers.hpp"
#include <iostream>

void SM(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
        double *Updates, unsigned int *Updates_index) {

  double C[Dim];
  double D[Dim];

  // For each update
  for (unsigned int l = 0; l < N_updates; l++) {

    // C = A^{-1} x U_l
    for (unsigned int i = 0; i < Dim; i++) {
      C[i] = 0;
      for (unsigned int j = 0; j < Dim; j++) {
        C[i] += Slater_inv[i * Dim + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (den < 1e-6) {
      std::cerr << "Breakdown condition triggered" << std::endl;
    }
    double iden = 1 / den;

    // D = v^T x A^{-1}
    for (unsigned int j = 0; j < Dim; j++) {
      D[j] = Slater_inv[(Updates_index[l] - 1) * Dim + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (unsigned int i = 0; i < Dim; i++) {
      for (unsigned int j = 0; j < Dim; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * Dim + j] -= update;
      }
    }
  }
}
