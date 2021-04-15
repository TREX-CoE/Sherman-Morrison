// SM-Standard.cpp
// Standard Sherman Morrison with multiple updates
#include "SM_Standard.hpp"
#include "SM_Helpers.hpp"

// Set common break-down threshold
double threshold() {
  double const threshold = 1e-3;
  return threshold;
}

// Naïve Sherman Morrison
void SM1(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index)
{

  std::cerr << "Called SM1 with " << N_updates << " updates" << std::endl;
  double C[Dim];
  double D[Dim];

  unsigned int l = 0;
  // For each update
  while (l < N_updates)
  {
    // C = A^{-1} x U_l
    for (unsigned int i = 0; i < Dim; i++)
    {
      C[i] = 0;
      for (unsigned int j = 0; j < Dim; j++)
      {
        C[i] += Slater_inv[i * Dim + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (fabs(den) < threshold())
    {
      std::cerr << "Breakdown condition triggered at " << Updates_index[l] << std::endl;
    }
    double iden = 1 / den;

    // D = v^T x A^{-1}
    for (unsigned int j = 0; j < Dim; j++)
    {
      D[j] = Slater_inv[(Updates_index[l] - 1) * Dim + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (unsigned int i = 0; i < Dim; i++)
    {
      for (unsigned int j = 0; j < Dim; j++)
      {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * Dim + j] -= update;
      }
    }

    l += 1;
  }
}

// Sherman Morrison, with J. Slagel splitting
// http://hdl.handle.net/10919/52966
void SM2(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index)
{

  std::cerr << "Called SM2 with updates " << N_updates << std::endl;
  double C[Dim];
  double D[Dim];

  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;

  unsigned int l = 0;
  // For each update
  while (l < N_updates)
  {
    // C = A^{-1} x U_l
    for (unsigned int i = 0; i < Dim; i++)
    {
      C[i] = 0;
      for (unsigned int j = 0; j < Dim; j++)
      {
        C[i] += Slater_inv[i * Dim + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (fabs(den) < threshold())
    {
      std::cerr << "Breakdown condition triggered at " << Updates_index[l] << std::endl;

      // U_l = U_l / 2 (do the split)
      for (unsigned int j = 0; j < Dim; j++)
      {
        later_updates[later * Dim + j] = Updates[l * Dim + j] / 2.0;
        C[j] /= 2.0;
      }
      later_index[later] = Updates_index[l];
      later++;

      den = 1 + C[Updates_index[l] - 1];
    }
    double iden = 1 / den;

    // D = v^T x A^{-1}
    for (unsigned int j = 0; j < Dim; j++)
    {
      D[j] = Slater_inv[(Updates_index[l] - 1) * Dim + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (unsigned int i = 0; i < Dim; i++)
    {
      for (unsigned int j = 0; j < Dim; j++)
      {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * Dim + j] -= update;
      }
    }
    l += 1;
  }

  if (later > 0)
  {
    SM2(Slater_inv, Dim, later, later_updates, later_index);
  }
}

// Sherman Morrison, leaving zero denominators for later
void SM3(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index)
{

  std::cerr << "Called SM3 with updates " << N_updates << std::endl;
  double C[Dim];
  double D[Dim];

  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;

  unsigned int l = 0;
  // For each update
  while (l < N_updates)
  {
    // C = A^{-1} x U_l
    for (unsigned int i = 0; i < Dim; i++)
    {
      C[i] = 0;
      for (unsigned int j = 0; j < Dim; j++)
      {
        C[i] += Slater_inv[i * Dim + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (fabs(den) < threshold())
    {
      std::cerr << "Breakdown condition triggered at " << Updates_index[l] << std::endl;

      for (unsigned int j = 0; j < Dim; j++)
      {
        later_updates[later * Dim + j] = Updates[l * Dim + j];
      }
      later_index[later] = Updates_index[l];
      later++;
      l += 1;
      continue;
    }
    double iden = 1 / den;

    // D = v^T x A^{-1}
    for (unsigned int j = 0; j < Dim; j++)
    {
      D[j] = Slater_inv[(Updates_index[l] - 1) * Dim + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (unsigned int i = 0; i < Dim; i++)
    {
      for (unsigned int j = 0; j < Dim; j++)
      {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * Dim + j] -= update;
      }
    }
    l += 1;
  }

  // If all the updates have failed, exit early with an error
  if (later == N_updates) {
    std::cerr << "SM3 cannot invert this matrix" << std::endl;
    return;
  }
  if (later > 0)
  {
    SM3(Slater_inv, Dim, later, later_updates, later_index);
  }
}
