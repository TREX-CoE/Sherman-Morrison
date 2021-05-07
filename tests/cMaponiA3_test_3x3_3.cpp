// main.cpp
#include "SM_Maponi.hpp"
#include "SM_Helpers.hpp"

int main() {

    unsigned int M = 3; // Dimension of the Slater-matrix
    unsigned int i, j; // Indices for iterators

    // Declare, allocate all vectors and matrices and fill them with zeros
    unsigned int *Ar_index = new unsigned int[M];
    double *A = new double[M*M]; // The matrix to be inverted
    double *A0 = new double[M*M]; // A diagonal matrix with the digonal elements of A
    double *Ar = new double[M*M]; // The update matrix
    double *A0_inv = new double[M*M]; // The inverse

    // Fill with zeros
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            A0[i*M+j] = 0;
            Ar[i*M+j] = 0;
            A0_inv[i*M+j] = 0;
        }
    }

    // Initialize A with M=3 and fill acc. to Eq. (17) from paper
    A[0] = 1;    A[3] = 1;    A[6] = -1;
    A[1] = 1;    A[4] = 1;    A[7] = 0;
    A[2] = -1;   A[5] = 0;    A[8] = -1;

    showMatrix(A, M, "A");

    // Initialize the diagonal matrix A0,
    // the inverse of A0_inv of diagonal matrix A0_inv
    // and the update matrix Ar
    for (i = 0; i < M; i++) {
        A0[i*M + i] = A[i*M + i];
        A0_inv[i*M + i] = 1.0/A[i*M + i];
        Ar_index[i] = i+1; // ! First column needs to start with 1 !
        for (j = 0; j < M; j++) {
            Ar[i*M + j] = A[i*M + j] - A0[i*M + j];
        }
    }

    // Define pointers dim and n_updates to use in Sherman-Morrison(...) function call
    MaponiA3(A0_inv, M, M, Ar, Ar_index);
    showMatrix(A0_inv, M, "A0_inv");

    // Deallocate all vectors and matrices
    delete [] A, A0, A0_inv, Ar, Ar_index;

    return 0;
}
