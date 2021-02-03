// main.cpp
#include "SM-MaponiA3.hpp"
#include "Helpers.hpp"
#include <cstdlib>
#include <ctime>

int main() {

    srand((unsigned) time(0));
    unsigned int randRange = 1; // to get random integers in range [-randRange, randRange]
    unsigned int M = 3; // Dimension of the Slater-matrix
    unsigned int i, j; // Indices for iterators

    // Declare, allocate all vectors and matrices and fill them with zeros
    unsigned int *Ar_index = new unsigned int[M];
    int **A = new int*[M]; // The matrix to be inverted
    int **A0 = new int*[M]; // A diagonal matrix with the digonal elements of A
    int **Ar = new int*[M]; // The update matrix
    double **A0_inv = new double*[M]; // Inverse of A0
    for (i = 0; i < M; i++) { 
        A[i] = new int[M];
        A0[i] = new int[M];
        Ar[i] = new int[M];
        A0_inv[i] = new double[M];
    }
    // Fill with zeros
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            A0[i][j] = 0;
            Ar[i][j] = 0;
            A0_inv[i][j] = 0;
        }
    }

    // Initialize A with M=3 and fill acc. to Eq. (17) from paper
    A[0][0] = 1;    A[0][1] = 1;    A[0][2] = -1;
    A[1][0] = 1;    A[1][1] = 1;    A[1][2] = 0;
    A[2][0] = -1;   A[2][1] = 0;    A[2][2] = -1;
    // // Fill A with random numbers from [-randRange,randRange] 
    // // and check if A and A0 are invertable
    // do {
    //     for (i = 0; i < M; i++) {
    //         for (j = 0; j < M; j++) {
    //             A[i][j] = rand()%(2*randRange+1)-randRange;
    //         }
    //     }
    //     for (i = 0; i < M; i++) {
    //         A0[i][i] = A[i][i];
    //     }
    // } while (matDet(A, M) == 0 || matDet(A0, M) == 0);
    showMatrix(A, M, "A");

    // Initialize the diagonal matrix A0,
    // the inverse of A0_inv of diagonal matrix A0_inv
    // and the update matrix Ar
    for (i = 0; i < M; i++) {
        A0[i][i] = A[i][i];
        A0_inv[i][i] = 1.0/A[i][i];
        Ar_index[i] = i;
        for (j = 0; j < M; j++) {
            Ar[i][j] = A[i][j] - A0[i][j];
        }
    }

    // Define pointers dim and n_updates to use in Sherman-Morrison(...) function call
    unsigned int *dim = new unsigned int(M);
    unsigned int *n_updates = new unsigned int(M);
    Sherman_Morrison(A0, A0_inv, dim, n_updates, Ar, Ar_index);
    showMatrix(A0_inv, M, "A0_inv");
    
    // Deallocate all vectors and matrices 
    for (i = 0; i < M; i++) {
        delete [] A[i], A0[i], A0_inv[i], Ar[i];
    }
    delete [] A, A0, A0_inv, Ar, Ar_index;
    delete dim, n_updates;

    return 0;
}