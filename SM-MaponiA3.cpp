// SM-MaponiA3.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM-MaponiA3.hpp"
#include "Helpers.hpp"

void Sherman_Morrison(int **Slater0, double **Slater_inv, unsigned int *Dim, unsigned int *N_updates, int **Updates, unsigned int *Updates_index) {
    unsigned int k, l, lbar, i, j, tmp, M = *Dim;
    unsigned int *p = new unsigned int[M+1];
    double *breakdown = new double[M+1];
    double alpha, beta;

    for (i = 0; i < M+1; i++) {
        p[i] = i;
    }

    int **Id = new int*[M];
    for (i = 0; i < M; i++) Id[i] = new int[M];
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            if (i != j) Id[i][j] = 0;
            else Id[i][j] = 1;
        }
    }

    // Declare auxiliary solution matrix ylk
    double ***ylk = new double**[M];
    for (l = 0; l < M; l++) {
        ylk[l] = new double*[M+1];
        for (k = 0; k < M+1; k++) {
            ylk[l][k] = new double[M+1];
        }
    }
    // Initialize ylk with zeros
    for (l = 0; l < M; l++) {
        for (k = 0; k < M+1; k++) {
            for (i = 0; i < M+1; i++) {
                ylk[l][k][i] = 0;
            }
        }
    }

    // Calculate all the y0k in M^2 multiplications instead of M^3
    for (k = 1; k < M+1; k++) {
        for (i = 1; i < M+1; i++) {
            ylk[0][k][i] = Slater_inv[i-1][i-1] * Updates[i-1][k-1];
        }
    }

    // Calculate all the ylk from the y0k
    for (l = 1; l < M; l++) {
        for (j = l; j < M+1; j++) {
            breakdown[j] = abs( 1 + ylk[l-1][p[j]][p[j]] );
        }
        lbar = getMaxIndex(breakdown, M+1);
        for (i = 0; i < M; i++) {
            breakdown[i] = 0;
        }
        tmp = p[l];
        p[l] = p[lbar];
        p[lbar] = tmp;
        for (k = l+1; k < M+1; k++) {
            beta = 1 + ylk[l-1][p[l]][p[l]];
            if (beta == 0) {
                cout << "Break-down condition occured. Exiting..." << endl;
                exit;
            }
            for (i = 1; i < M+1; i++) {
                alpha = ylk[l-1][p[k]][p[l]] / beta;
                ylk[l][p[k]][i] = ylk[l-1][p[k]][i] - alpha * ylk[l-1][p[l]][i];
            }
        }
    }

    // Construct A-inverse from A0-inverse and the ylk
    double **U;
    double **Al = new double*[M];
    for (i = 0; i < M; i++) Al[i] = new double[M];
  
    for (l = 0; l < M; l++) {
        k = l+1;
        U = outProd(ylk[l][p[k]], Id[p[k]-1], M);
        beta = 1 + ylk[l][p[k]][p[k]];
        for (i = 0; i < M; i++) {
            for (j = 0; j < M; j++) {
                Al[i][j] = Id[i][j] - U[i][j] / beta;
            }
        }
        showMatrix(Slater_inv, M, "Slater_inv");
        Slater_inv = matMul(Al, Slater_inv, M);
        showMatrix(Slater_inv, M, "Slater_inv");
    }

    delete [] p, breakdown;

    for (i = 0; i < M; i++) {
        delete [] Id[i];
        delete [] U[i];
        delete [] Al[i];
    }

    for (l = 0; l < M; l++) {
        for (k = 0; k < M+1; k++) {
            delete [] ylk[l][k];
        }
        delete [] ylk[l];
    }
}