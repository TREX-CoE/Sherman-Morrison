// SM-MaponiA3_f.cpp
// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007
#include "SM_MaponiA3_f.hpp"
#include "Helpers.hpp"

void MaponiA3(int **linSlater0, double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates, int **linUpdates, unsigned int *Updates_index) {

    // Define new 2D arrays and copy the elements of the 
    // linear passed Fortran arrays. This block needs to
    // be replaced with some casting mechanism to avoid
    // copying of arrays.
    int **Slater0 = new int*[*Dim];
    int **Updates = new int*[*Dim];
    double **Slater_inv = new double*[*Dim];
    for (int i = 0; i < *Dim; i++) {
        Slater0[i] = new int[*Dim];
        Updates[i] = new int[*Dim];
        Slater_inv[i] = new double[*Dim];
    }
    for (unsigned int i = 0; i < *Dim; i++) {
        for (unsigned int j = 0; j < *Dim; j++) {
            Slater0[i][j] = linSlater0[0][i+*Dim*j];
            Slater_inv[i][j] = linSlater_inv[0][i+*Dim*j];
            Updates[i][j] = linUpdates[0][i+*Dim*j];
        }
    }

    // Possible casting candidates
    // int (*Slater0)[*Dim] = (int(*)[*Dim])linSlater0[0];
    // double (*Slater_inv)[*Dim] = (double(*)[*Dim])linSlater_inv[0];
    // int (*Updates)[*Dim] = (int(*)[*Dim])linUpdates[0];
    ////////////////////////////////////////////////////////////////////////

    unsigned int k, l, lbar, i, j, tmp, M = *Dim;
    unsigned int *p = new unsigned int[M+1];
    unsigned int **Id = new unsigned int*[M];
    double alpha, beta;
    double **U, *breakdown = new double[M+1];
    double **Al = new double*[M];
    p[0] = 0;
    for (i = 0; i < M; i++) {
        p[i+1] = i + 1;
        Id[i] = new unsigned int[M];
        Al[i] = new double[M];
    }

    // Declare auxiliary solution matrix ylk
    double ***ylk = new double**[M];
    for (l = 0; l < M; l++) {
        ylk[l] = new double*[M+1];
        for (k = 0; k < M+1; k++) {
            ylk[l][k] = new double[M+1];
        }
    }

    // Initialize identity matrix
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            if (i != j) Id[i][j] = 0;
            else Id[i][j] = 1;
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
    // Keep the memory location of the passed array 'Slater_inv' before 'Slater_inv'
    // gets reassigned by 'matMul(...)' in the next line, by creating a new
    // pointer 'copy' that points to whereever 'Slater_inv' points to now.
    // double **copy  = Slater_inv;

    for (l = 0; l < M; l++) {
        k = l+1;
        U = outProd(ylk[l][p[k]], Id[p[k]-1], M);
        beta = 1 + ylk[l][p[k]][p[k]];
        for (i = 0; i < M; i++) {
            for (j = 0; j < M; j++) {
                Al[i][j] = Id[i][j] - U[i][j] / beta;
            }
        }
        Slater_inv = matMul(Al, Slater_inv, M);
    }

    // // Assign the new values of 'Slater_inv' to the old values in 'copy[][]'
    // for (i = 0; i < M; i++) {
    //     for (j = 0; j < M; j++) {
    //         copy[i][j] = Slater_inv[i][j];
    //     }
    // }

    // Assign the new values of 'Slater_inv' to the old values in 'copy[][]'
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            linSlater_inv[0][i+*Dim*j] = Slater_inv[i][j];
        }
    }

    for (l = 0; l < M; l++) {
        for (k = 0; k < M+1; k++) {
            delete [] ylk[l][k];
        }
        delete [] ylk[l], Id[l], U[l], Al[l], Slater_inv[l];
    }
    delete [] p, breakdown;
}