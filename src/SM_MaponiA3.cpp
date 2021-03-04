// SM-MaponiA3_f.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM_MaponiA3.hpp"
#include "Helpers.hpp"

void MaponiA3(double *Slater_inv, unsigned int Dim,
              unsigned int N_updates, double *Updates,
              unsigned int *Updates_index) {

  /*
  DECLARE AND INITIALISE ARRAYS
  */

  unsigned int k, l, lbar, i, j, tmp, component, max, breakdown;
  unsigned int *p = new unsigned int[N_updates + 1] {0};
  double alpha, beta;
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
  
  /*
  START THE ALGORITHM
  */

  // Calculate the y0k
  for (k = 1; k < N_updates + 1; k++) {
    for (i = 1; i < Dim + 1; i++) {
      for (j = 1; j < Dim + 1; j++) {
        ylk[0][k][i] += Slater_inv[(i-1)*Dim + (j-1)]
                      * Updates[(k-1)*Dim + (j-1)];
      }
    }
  }

  // // Compute A_1_inv and A_1_inv * Slater_inv
  // double *lastIntermRes = Slater_inv;
  // double *Res = new double[Dim*Dim]{0};
  // component = Updates_index[0];
  // beta = 1 + ylk[0][1][component];
  // for (i = 0; i < Dim; i++) {
  //   for (j = 0; j < Dim; j++) {
  //     Al[i*Dim + j] = (i == j) - (j   == component-1) * ylk[0][1][i + 1] / beta;
  //   }
  // }
  // matMul2(Al, Slater_inv, Res, Dim);
  // Slater_inv = Res;
  // Res = lastIntermRes;
  // lastIntermRes = Slater_inv;
  // showMatrix(Slater_inv, Dim, "Slater_inv");

  // Calculate all the ylk from the y0k
  for (l = 0; l < N_updates; l++) {

    // Select update with largest break-down val
    lbar = l+1; max = 0;
    for (j = l+1; j < N_updates + 1; j++) {
      component = Updates_index[p[j] - 1];
      breakdown = abs(1 + ylk[l+1 - 1][p[j]][component]);
      if (breakdown > max) {
          max = breakdown;
          lbar = j;
        }
    } tmp = p[l+1]; p[l+1] = p[lbar]; p[lbar] = tmp;

    component = Updates_index[p[l+1] - 1];
    beta = 1 + ylk[l][p[l+1]][component];
    if (beta == 0) {
      cout << "Break-down occured. Exiting..." << endl;
      exit(1);
    }
    for (k = l+2; k < N_updates + 1; k++) {
      alpha = ylk[l][p[k]][component] / beta;
      cout << "( l, k, p[k], component ) = (" << l << ", " << k << ", " << p[k] << ", " << component << ")" << endl;
      for (i = 1; i < Dim + 1; i++) {
        cout << "ylk[" << l << "][p[" << k << "]][" << i << "] = ylk[" << l - 1 << "][p[" << k << "]][" << i << "] - alpha * ylk[" << l - 1 << "][p[" << l << "]][" << i << "]" << endl;
        ylk[l+1][p[k]][i] = ylk[l][p[k]][i]
                        - alpha * ylk[l][p[l+1]][i];
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
  void MaponiA3_f(double **linSlater_inv, unsigned int *Dim,
                  unsigned int *N_updates, double **linUpdates,
                  unsigned int **Updates_index) {
    MaponiA3(*linSlater_inv, *Dim,
             *N_updates, *linUpdates,
             *Updates_index);
  }
}
