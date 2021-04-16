// SM-MaponiA3_f.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM_MaponiA3S.hpp"
#include "SM_Helpers.hpp"

// #define DEBUG

void MaponiA3S(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
               double *Updates, unsigned int *Updates_index) {
  /*
  DECLARE AND INITIALISE ARRAYS
  */

  unsigned int k, l, i, j, component;
  unsigned int *p = new unsigned int[N_updates + 1]{0};
  double alpha, beta;
  double *Al = new double[Dim * Dim];
  double *next = new double[Dim * Dim]{0};
  double *last = Slater_inv, *tmp;

  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;

  // Populate update-order vector
  for (i = 0; i < N_updates; i++) {
    p[i + 1] = i + 1;
  }

  // Declare auxiliary solution matrix ylk[N_updates][N_updates+1][Dim]
  double ***ylk = new double **[N_updates];
  for (l = 0; l < N_updates; l++) {
    ylk[l] = new double *[N_updates + 1];
    for (k = 0; k < N_updates + 1; k++) {
      ylk[l][k] = new double[Dim + 1]{0};
    }
  }

  /*
  START ALGORITHM
  */

  // Calculate the {y_{0,k}}
  for (k = 1; k < N_updates + 1; k++) {
#ifdef DEBUG
    std::cout << "Compute y0k: " << std::endl;
    std::cout << "ylk[0][" << k << "][:]";
    std::cout << std::endl;
#endif
    for (i = 1; i < Dim + 1; i++) {
      for (j = 1; j < Dim + 1; j++) {
        ylk[0][k][i] += Slater_inv[(i - 1) * Dim + (j - 1)] *
                        Updates[(k - 1) * Dim + (j - 1)];
      }
    }
#ifdef DEBUG
    showVector(ylk[0][k], Dim, "");
#endif
  }

  // Calculate the {y_{l,k}} from the {y_{0,k}}
  for (l = 0; l < N_updates; l++) {
#ifdef DEBUG
    std::cout << "In outer compute-ylk-loop: l = " << l << std::endl;
    std::cout << std::endl;
#endif

    // For given l select intermediate update with largest break-down val
    selectLargestDenominator(l, N_updates, Updates_index, p, ylk);

    // Select component and comp. bd-condition.
    component = Updates_index[p[l + 1] - 1];
    beta = 1 + ylk[l][p[l + 1]][component];
#ifdef DEBUG
    std::cout << "p[l+1] = " << p[l + 1] << std::endl;
    std::cout << "component = " << component << std::endl;
    std::cout << "beta = 1 + ylk[" << l << "][" << p[l + 1] << "][" << component
              << "] = " << beta << std::endl;
    std::cout << std::endl;
#endif
    if (fabs(beta) < threshold()) {
      std::cerr << "Breakdown condition triggered at " << component
                << std::endl;
      
      for (unsigned int i = 0; i < Dim; i++) {
        later_updates[later * Dim + i] = Updates[l * Dim + i] * 0.5;
        ylk[l][p[l + 1]][i] *= 0.5;
      }
      later_index[later] = Updates_index[p[l + 1]];
      later++;
      beta = 1 + ylk[l][p[l + 1]][component];
    }
    double ibeta = 1.0 / beta;

// Compute intermediate update to Slater_inv
#ifdef DEBUG
    std::cout << "Compute intermediate update to Slater_inv" << std::endl;
    std::cout << "component = " << component << std::endl;
    std::cout << "beta = 1 + ylk[" << l << "][" << p[l + 1] << "][" << component
              << "]" << std::endl;
    std::cout << "ylk[l][p[k]][:] = ylk[" << l << "][" << p[l + 1] << "][:]"
              << std::endl;
    std::cout << std::endl;
#endif
    for (i = 0; i < Dim; i++) {
      for (j = 0; j < Dim; j++) {
        Al[i * Dim + j] =
            (i == j) - (j == component - 1) * ylk[l][p[l + 1]][i + 1] * ibeta;
      }
    }
    matMul(Al, last, next, Dim); tmp = next; next = last; last = tmp;
#ifdef DEBUG
    showMatrix(last, Dim, "last");
#endif

    // For given l != 0 compute the next {y_{l,k}}
    for (k = l + 2; k < N_updates + 1; k++) {
      alpha = ylk[l][p[k]][component] * ibeta;
#ifdef DEBUG
      std::cout << "Inside k-loop: k = " << k << std::endl;
      std::cout << "ylk[" << l + 1 << "][" << p[k] << "][:]" << std::endl;
      std::cout << std::endl;
#endif
      for (i = 1; i < Dim + 1; i++) {
        ylk[l + 1][p[k]][i] = ylk[l][p[k]][i] - alpha * ylk[l][p[l + 1]][i];
      }
    }
  }
  memcpy(Slater_inv, last, Dim * Dim * sizeof(double));

  if (later > 0) {
    MaponiA3S(Slater_inv, Dim, later, later_updates, later_index);
  }

  /*
  CLEANUP MEMORY
  */

  for (l = 0; l < N_updates; l++) {
    for (k = 0; k < N_updates + 1; k++) {
      delete[] ylk[l][k];
    }
    delete[] ylk[l];
  }
  delete[] Al, next, p;
}

extern "C" {
  void MaponiA3S_f(double **linSlater_inv, unsigned int *Dim,
                  unsigned int *N_updates, double **linUpdates,
                  unsigned int **Updates_index) {
    MaponiA3S(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index);
  }
}
