// SM-MaponiA3_f.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM_MaponiA3.hpp"
#include "Helpers.hpp"

void MaponiA3(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
              double *Updates, unsigned int *Updates_index) {

  unsigned int k, l, lbar, i, j, tmp, component;
  unsigned int *p = new unsigned int[N_updates + 1];
  double alpha, beta;
  double *breakdown = new double[N_updates + 1];
  double *Al = new double[Dim * Dim];
  p[0] = 0;
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
        ylk[0][k][i] += Slater_inv[(i-1) + (j-1)*Dim]
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
    // Reset breakdown back to 0 for next round to avoid case where its
    // first element is always the largest
    for (i = 0; i < N_updates + 1; i++) {
      breakdown[i] = 0;
    }
    tmp = p[l];
    p[l] = p[lbar];
    p[lbar] = tmp;
    component = Updates_index[p[l] - 1];
    beta = 1 + ylk[l - 1][p[l]][component];
    if (beta == 0) {
      cout << "Break-down occured. Exiting..." << endl;
      exit;
    }
    for (k = l + 1; k < N_updates + 1; k++) {
      alpha = ylk[l - 1][p[k]][component] / beta;
      for (i = 1; i < Dim + 1; i++) {
        ylk[l][p[k]][i] = ylk[l - 1][p[k]][i] - alpha * ylk[l - 1][p[l]][i];
      }
    }
  }

  // Keep the memory location of the passed array 'Slater_inv' before
  // 'Slater_inv' gets reassigned by 'matMul(...)' in the next line, by creating
  // a new pointer 'copy' that points to whereever 'Slater_inv' points to now.
  double *copy = Slater_inv;

  Slater_inv = transpose(Slater_inv, Dim);

  // Construct A-inverse from A0-inverse and the ylk
  for (l = 0; l < N_updates; l++) {       // l    = 0, 1
    k = l + 1;                            // k    = 1, 2
    component = Updates_index[p[k] - 1];  // comp = 2, 4
    beta = 1 + ylk[l][p[k]][component];
    for (i = 0; i < Dim; i++) {
      for (j = 0; j < Dim; j++) {
        Al[i * Dim + j] = (i == j) - (j == component - 1) * ylk[l][p[k]][i + 1] / beta;
      }
    }
    Slater_inv = matMul(Al, Slater_inv, Dim);
  }

  Slater_inv = transpose(Slater_inv, Dim);

  // Assign the new values of 'Slater_inv' to the old values in 'copy[][]'
  for (i = 0; i < Dim; i++) {
    for (j = 0; j < Dim; j++) {
      copy[i * Dim + j] = Slater_inv[i * Dim + j];
    }
  }

  for (l = 0; l < N_updates; l++) {
    for (k = 0; k < N_updates + 1; k++) {
      delete[] ylk[l][k];
    }
    delete[] ylk[l];
  }
  delete[] Al;
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
