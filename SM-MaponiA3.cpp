// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007

#include <iostream>
using namespace std;

int M = 3;

int main() {

    uint i, j, k, l, lbar, tmp;
    // Declare and allocate all vectors and matrices
    int* p = new int[M+1];
    int* breakdown = new int[M];

    int** A = new int*[M];
    int** A0 = new int*[M];
    int** Ar = new int*[M];
    int** Id = new int*[M];
    double** A0inv = new double*[M];
    double** Ainv = new double*[M];
    for (i = 0; i < M; i++) { 
        A[i] = new int[M];
        Ainv[i] = new double[M];
        A0[i] = new int[M];
        A0inv[i] = new double[M];
        Ar[i] = new int[M];
        Id[i] = new int[M];
    }

    double*** ylk = new double**[M];
    for (l = 0; l < M; l++) {
        ylk[l] = new double*[M+1];
        for (k = 0; k < M+1; k++) {
            ylk[l][k] = new double[M];
        }
    }

    // Initialize all matrices with zeros
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            A[i][j] = 0;
            Ainv[i][j] = 0;
            A0[i][j] = 0;
            A0inv[i][j] = 0;
            Ar[i][j] = 0;
            Id[i][j] = 0;
        }
    }

    // Initialize all vectors and matrices
    A[0][0] = 1;    A[0][1] = 1;    A[0][2] = -1;
    A[1][0] = 1;    A[1][1] = 1;    A[1][2] = 0;
    A[2][0] = -1;   A[2][1] = 0;    A[2][2] = -1;

    // Define identity matrix, A0, A0inv and p
    for (i = 0; i < M; i++) {
        Id[i][i] = 1;
        A0[i][i] = A[i][i];
        A0inv[i][i] = 1.0/A[i][i];
        p[i] = i;
    }

    // Init Ar
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            Ar[i][j] = A[i][j] - A0[i][j];
        }
    }

    // Calculate all the y0k in M^2 multiplications instead of M^3
    for (k = 1; i < M+1; k++) {
        for (i = 0; i < M; i++) {
            ylk[0][k][i] = A0inv[i][i] * Ar[i][k];
        }
    }

    // Calculate all the ylk from the y0k
    for (l = 1; l < M; l++) {
        for (j = l; j < M; j++) {
            breakdown[j] = abs( 1 + ylk[l][p[j]][p[j]-1] );
        }
        lbar = findMaxElmnt(breakdown); // NOT IMPLEMENTED YET
        for (i = 0; i < M; i++) {
            breakdown[i] = 0;
        }
        tmp = p
    }



    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            cout << "A["<<i<<"]["<<j<<"] = " << A0inv[i][j] << " ";
        }
        cout << endl;
    }







    // Deallocate all vectors and matrices 
    for (i = 0; i < M; i++) {
        delete [] A[i];
        delete [] Ainv[i];
        delete [] A0[i];
        delete [] A0inv[i];
        delete [] Ar[i];
        delete [] Id[i];
    }

    for (l = 0; l < M; l++) {
        for (k = 0; k < M+1; k++) {
            delete [] ylk[l][k];
        }
        delete [] ylk[l];
    }

    delete [] p, breakdown, A, Ainv, A0, A0inv, Ar, Id, ylk;

    return 0;
}
