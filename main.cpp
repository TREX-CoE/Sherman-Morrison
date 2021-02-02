// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007

#include <cstdlib>
#include <ctime>
#include "SM-MaponiA3.hpp"
using namespace std;

int main() {

    srand((unsigned) time(0));
    uint randRange = 1; // to get random integers in range [-randRange, randRange]
    uint M = 3; // Dimension of the Slater-matrix
    uint i, j; // Indices for iterators

    // Declare and allocate all vectors and matrices
    uint *Ar_index = new uint[M];
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

    // Initialize all matrices with zeros
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            A0[i][j] = 0;
            Ar[i][j] = 0;
            A0_inv[i][j] = 0;
        }
    }

    // Initialize A with M=3 and Eq. (17) from paper
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

    // Init the update matrix Ar, A0_inv and Ar_index
    for (i = 0; i < M; i++) {
        A0_inv[i][i] = 1.0/A[i][i];
        Ar_index[i] = i;
        for (j = 0; j < M; j++) {
            Ar[i][j] = A[i][j] - A0[i][j];
        }
    }

    uint *dim = new uint(M);
    uint *n_updates = new uint(M);
    Sherman_Morrison(A0, A0_inv, dim, n_updates, Ar, Ar_index);
    
    showMatrix(A0_inv, M, "A0_inv");
    
    // Deallocate all vectors and matrices 
    for (i = 0; i < M; i++) {
        delete [] A[i];
        delete [] A0[i];
        delete [] A0_inv[i];
        delete [] Ar[i];
    }
    delete [] A, A0, A0_inv, Ar, Ar_index;
    delete dim, n_updates;

    return 0;
}
