// SM-MaponiA3_f.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM_MaponiA3.hpp"
#include "Helpers.hpp"

void MaponiA3(double *Slater_inv, unsigned int Dim,
              unsigned int N_updates, double *Updates,
              unsigned int *Updates_index) {

  unsigned int k, l, lbar, i, j, tmp, component;
  unsigned int *p = new unsigned int[N_updates + 1] {0};
  double alpha, beta;
  double *breakdown = new double[N_updates + 1] {0};
  double *Al = new double[Dim * Dim];

  // Populate update-order vector
  for (i = 0; i < N_updates; i++) {
    p[i + 1] = i + 1;
  }

  // Declare auxiliary solution matrix ylk
  double ***ylk = new double **[N_updates];
  for (l = 0; l < N_updates; l++) {
    ylk[l] = new double *[N_updates + 1];
    for (k = 0; k < N_updates + 1; k++) {
      ylk[l][k] = new double[Dim + 1] {0};
    }
  }

  // Calculate the y0k
  for (k = 1; k < N_updates + 1; k++) {
    for (i = 1; i < Dim + 1; i++) {
      for (j = 1; j < Dim + 1; j++) {
        ylk[0][k][i] += Slater_inv[(i-1)*Dim + (j-1)]
                      * Updates[(k-1)*Dim + (j-1)];
      }
    }
  }

  // Calculate all the ylk from the y0k
  for (l = 1; l < N_updates; l++) {
    for (j = l; j < N_updates + 1; j++) {
      component = Updates_index[p[j] - 1];
      breakdown[j] = abs(1 + ylk[l - 1][p[j]][component]);
    }
    lbar = getMaxIndex(breakdown, N_updates + 1);
    // Reset breakdown back to 0 for next round to avoid case where
    // its first element is always the largest
    for (i = 0; i < N_updates + 1; i++) {
      breakdown[i] = 0;
    }
    tmp = p[l];
    p[l] = p[lbar];
    p[lbar] = tmp;
    component = Updates_index[p[l] - 1];
    beta = 1 + ylk[l - 1][p[l]][component];
    if (fabs(beta) < 1e-6) {
      cout << "Break-down occured. Exiting..." << endl;
      exit(1);
    }
    for (k = l + 1; k < N_updates + 1; k++) {
      alpha = ylk[l - 1][p[k]][component] / beta;
      cout << "( l, k, p[k], component ) = (" << l << ", " << k << ", " << p[k] << ", " << component << ")" << endl;
      for (i = 1; i < Dim + 1; i++) {
        cout << "ylk[" << l << "][p[" << k << "]][" << i << "] = ylk[" << l - 1 << "][p[" << k << "]][" << i << "] - alpha * ylk[" << l - 1 << "][p[" << l << "]][" << i << "]" << endl;
        ylk[l][p[k]][i] = ylk[l - 1][p[k]][i]
                        - alpha * ylk[l - 1][p[l]][i];
      }
    }

  }
  
  // Construct A-inverse from A0-inverse and the ylk
  double *last = Slater_inv;
  double *next = new double[Dim*Dim] {0};
  for (l = 0; l < N_updates; l++) {
    k = l + 1;
    component = Updates_index[p[k] - 1];
    beta = 1 + ylk[l][p[k]][component];
    cout << "( l, k, p[k], component ) = (" << l << ", " << k << ", " << p[k] << ", " << component << ")" << endl;
    cout << "ylk[" << l << "][" << p[k] << "][i + 1]" << endl;
    for (i = 0; i < Dim; i++) {
      for (j = 0; j < Dim; j++) {
        Al[i*Dim + j] = (i == j) - (j == component-1)
                      * ylk[l][p[k]][i + 1] / beta;
      }
    }
    matMul2(Al, last, next, Dim);
    double *tmp = next;
    next = last;
    last = tmp;
  }
  memcpy(Slater_inv, last, Dim*Dim*sizeof(double));
 
  // Free memory
  for (l = 0; l < N_updates; l++) {
    for (k = 0; k < N_updates + 1; k++) {
      delete[] ylk[l][k];
    }
    delete[] ylk[l];
  }
  delete[] Al, next;
  delete[] p, breakdown;
}

extern "C" {
  void MaponiA3_f(double **linSlater_inv, unsigned int *Dim,
                  unsigned int *N_updates, double **linUpdates,
                  unsigned int **Updates_index) {
    MaponiA3(*linSlater_inv, *Dim,
             *N_updates, *linUpdates,
             *Updates_index);
  }
}
