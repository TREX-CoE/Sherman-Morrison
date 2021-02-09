// SM-MaponiA3_f.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM_MaponiA3.hpp"
#include "Helpers.hpp"

void MaponiA3(double *Slater0, double *Slater_inv, unsigned int M,
              unsigned int N_updates, double *Updates,
              unsigned int *Updates_index) {

  unsigned int k, l, lbar, i, j, tmp = M;
  unsigned int *p = new unsigned int[M + 1];
  double alpha, beta;
  double *breakdown = new double[M + 1];
  double *Al = new double[M * M];
  p[0] = 0;
  for (i = 0; i < M; i++) {
    p[i + 1] = i + 1;
  }

  // Declare auxiliary solution matrix ylk
  double ***ylk = new double **[M];
  for (l = 0; l < M; l++) {
    ylk[l] = new double *[M + 1];
    for (k = 0; k < M + 1; k++) {
      ylk[l][k] = new double[M + 1];
    }
  }

  // Initialize ylk with zeros
  for (l = 0; l < M; l++) {
    for (k = 0; k < M + 1; k++) {
      for (i = 0; i < M + 1; i++) {
        ylk[l][k][i] = 0;
      }
    }
  }

  // Calculate all the y0k in M^2 multiplications instead of M^3
  for (k = 1; k < M + 1; k++) {
    for (i = 1; i < M + 1; i++) {
      ylk[0][k][i] =
          Slater_inv[(i - 1) * M + (i - 1)] * Updates[(i - 1) * M + (k - 1)];
    }
  }

  // Calculate all the ylk from the y0k
  for (l = 1; l < M; l++) {
    for (j = l; j < M + 1; j++) {
      breakdown[j] = abs(1 + ylk[l - 1][p[j]][p[j]]);
    }
    lbar = getMaxIndex(breakdown, M + 1);
    for (i = 0; i < M; i++) {
      breakdown[i] = 0;
    }
    tmp = p[l];
    p[l] = p[lbar];
    p[lbar] = tmp;
    for (k = l + 1; k < M + 1; k++) {
      beta = 1 + ylk[l - 1][p[l]][p[l]];
      if (beta == 0) {
        cout << "Break-down condition occured. Exiting..." << endl;
        exit;
      }
      for (i = 1; i < M + 1; i++) {
        alpha = ylk[l - 1][p[k]][p[l]] / beta;
        ylk[l][p[k]][i] = ylk[l - 1][p[k]][i] - alpha * ylk[l - 1][p[l]][i];
      }
    }
  }

  // Keep the memory location of the passed array 'Slater_inv' before
  // 'Slater_inv' gets reassigned by 'matMul(...)' in the next line, by creating
  // a new pointer 'copy' that points to whereever 'Slater_inv' points to now.
  double *copy = Slater_inv;

  double *U = new double[M * M];
  // Construct A-inverse from A0-inverse and the ylk
  for (l = 0; l < M; l++) {
    k = l + 1;
    for (unsigned int i = 0; i < M; i++) {
      for (unsigned int j = 0; j < M; j++) {
        U[i * M + j] = ylk[l][p[k]][i + 1] * ((p[k] - 1) == j);
      }
    }

    beta = 1 + ylk[l][p[k]][p[k]];
    for (i = 0; i < M; i++) {
      for (j = 0; j < M; j++) {
        Al[i * M + j] = (i == j) - U[i * M + j] / beta;
      }
    }
    Slater_inv = matMul(Al, Slater_inv, M);
  }

  // Assign the new values of 'Slater_inv' to the old values in 'copy[][]'
  for (i = 0; i < M; i++) {
    for (j = 0; j < M; j++) {
      copy[i * M + j] = Slater_inv[i * M + j];
    }
  }

  for (l = 0; l < M; l++) {
    for (k = 0; k < M + 1; k++) {
      delete[] ylk[l][k];
    }
    delete[] ylk[l];
  }
  delete[] Al, U;
  delete[] p, breakdown;
}

extern "C" {
void MaponiA3_f(double **linSlater0, double **linSlater_inv, unsigned int *Dim,
                unsigned int *N_updates, double **linUpdates,
                unsigned int **Updates_index) {
  MaponiA3(*linSlater0, *linSlater_inv, *Dim, *N_updates, *linUpdates,
           *Updates_index);
}
}
