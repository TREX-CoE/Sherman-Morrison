// Algorithm 3 from P. Maponi,
// p. 283, doi:10.1016/j.laa.2006.07.007

#include <iostream>
#include <string>
using namespace std;

uint getMaxIndex(double* arr, uint size);
void showScalar(uint scalar, string name);
void showScalar(int scalar, string name);
void showScalar(double scalar, string name);
void showVector(uint* vector, uint size, string name);
void showVector(int* vector, uint size, string name);
void showVector(double* vector, uint size, string name);
void showMatrix(uint** matrix, uint sizeR, uint sizeC, string name);
void showMatrix(int** matrix, uint sizeR, uint sizeC, string name);
void showMatrix(double** matrix, uint sizeR, uint sizeC, string name);
void showMatrixT(uint** matrix, uint sizeR, uint sizeC, string name);
void showMatrixT(int** matrix, uint sizeR, uint sizeC, string name);
void showMatrixT(double** matrix, uint sizeR, uint sizeC, string name);
int** matMul(int** A, int** B, uint size);
double** matMul(double** A, double** B, uint size);
int** outProd(int* vec1, int* vec2, uint size);
double** outProd(double* vec1, int* vec2, uint size);

int main() {

    uint M = 3;
    uint i, j, k, l, lbar, tmp;
    double alpha, beta;

    // Declare and allocate all vectors and matrices
    uint* p = new uint[M+1];
    double* breakdown = new double[M+1];

    int** A = new int*[M];
    int** A0 = new int*[M];
    int** Ar = new int*[M];
    int** Id = new int*[M];
    double** A0inv = new double*[M];
    double** Ainv; //= new double*[M];
    for (i = 0; i < M; i++) { 
        A[i] = new int[M];
        //Ainv[i] = new double[M];
        A0[i] = new int[M];
        A0inv[i] = new double[M];
        Ar[i] = new int[M];
        Id[i] = new int[M];
    }

    double*** ylk = new double**[M];
    for (l = 0; l < M; l++) {
        ylk[l] = new double*[M+1];
        for (k = 0; k < M+1; k++) {
            ylk[l][k] = new double[M+1];
        }
    }

    // Initialize all matrices with zeros
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            A[i][j] = 0;
            //Ainv[i][j] = 0;
            A0[i][j] = 0;
            A0inv[i][j] = 0;
            Ar[i][j] = 0;
            Id[i][j] = 0;
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

    // Initialize A
    A[0][0] = 1;    A[0][1] = 1;    A[0][2] = -1;
    A[1][0] = 1;    A[1][1] = 1;    A[1][2] = 0;
    A[2][0] = -1;   A[2][1] = 0;    A[2][2] = -1;

    // Define identity matrix, A0, A0inv and p
    p[0] = 0;
    for (i = 0; i < M; i++) {
        Id[i][i] = 1;
        A0[i][i] = A[i][i];
        A0inv[i][i] = 1.0/A[i][i];
        p[i+1] = i+1;
    }

    // Init Ar
    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            Ar[i][j] = A[i][j] - A0[i][j];
        }
    }

    showMatrix(A, M, M, "A");
    // showMatrix(A0, M, M, "A0");
    // showMatrix(Ar, M, M, "Ar");
    // showMatrix(A0inv, M, M, "A0inv");

    // Calculate all the y0k in M^2 multiplications instead of M^3
    for (k = 1; k < M+1; k++) {
        for (i = 1; i < M+1; i++) {
            ylk[0][k][i] = A0inv[i-1][i-1] * Ar[i-1][k-1];
        }
        // showVector(ylk[0][k], M+1, "y0k");
    }
    showMatrixT(ylk[0], M+1, M+1, "y0k");

    // Calculate all the ylk from the y0k
    // showVector(p, M+1, "p");
    for (l = 1; l < M; l++) {
        for (j = l; j < M+1; j++) {
            breakdown[j] = abs( 1 + ylk[l-1][p[j]][p[j]] );
        }
        // showVector(breakdown, M+1, "break-down vector");
        lbar = getMaxIndex(breakdown, M+1);
        // showScalar(lbar, "lbar");
        for (i = 0; i < M; i++) {
            breakdown[i] = 0;
        }
        tmp = p[l];
        p[l] = p[lbar];
        p[lbar] = tmp;
        // showVector(p, M+1, "p");
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
            // showVector(ylk[l][p[k]], M+1, "ylk");
        }
    }
    showMatrixT(ylk[1], M+1, M+1, "y1k");
    showMatrixT(ylk[2], M+1, M+1, "y2k");
    // EVERYTHING WORKS UPTO HERE

    // Construct A-inverse from A0-inverse and the ylk
    double** U;
    double** Ap = new double*[M];
    for (i = 0; i < M; i++) Ap[i] = new double[M];
    Ainv = A0inv;
    for (l = 0; l < M; l++) {
        k = l+1;
        U = outProd(ylk[l][p[k]], Id[p[k]-1], M);
        beta = 1 + ylk[l][p[k]][p[k]];
        for (i = 0; i < M; i++) {
            for (j = 0; j < M; j++) {
                Ap[i][j] = Id[i][j] - U[i][j] / beta;
            }
        }
        Ainv = matMul(Ap, Ainv, M);
    }
    showMatrixT(Ainv, M, M, "Ainv");
    
    // Deallocate all vectors and matrices 
    for (i = 0; i < M; i++) {
        delete [] A[i];
        delete [] A0[i];
        delete [] A0inv[i];
        delete [] Ar[i];
        delete [] Id[i];
        delete [] U[i];
        delete [] Ap[i];
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

uint getMaxIndex(double* arr, uint size) { 
    uint i;     
    uint max = arr[0]; 
    uint maxi = 0;
    for (i = 1; i < size; i++) {
        if (arr[i] > max) {
            max = arr[i];
            maxi = i;
        }
    }
    return maxi; 
}

void showScalar(uint scalar, string name) {
    cout << name << " = " << scalar << endl << endl;
}
void showScalar(int scalar, string name) {
    cout << name << " = " << scalar << endl << endl;
}
void showScalar(double scalar, string name) {
    cout << name << " = " << scalar << endl << endl;
}

void showVector(uint* vector, uint size, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < size; i++) {
        cout << "[ " << vector[i] << " ]" << endl;
    }
    cout << endl;
}
void showVector(int* vector, uint size, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < size; i++) {
        cout << "[ " << vector[i] << " ]" << endl;
    }
    cout << endl;
}
void showVector(double* vector, uint size, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < size; i++) {
        cout << "[ " << vector[i] << " ]" << endl;
    }
    cout << endl;
}

void showMatrix(uint** matrix, uint sizeR, uint sizeC, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < sizeR; i++) {
        cout << "[ ";
        for (uint j = 0; j < sizeC; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}
void showMatrix(int** matrix, uint sizeR, uint sizeC, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < sizeR; i++) {
        cout << "[ ";
        for (uint j = 0; j < sizeC; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}

void showMatrix(double** matrix, uint sizeR, uint sizeC, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < sizeR; i++) {
        cout << "[ ";
        for (uint j = 0; j < sizeC; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}
void showMatrixT(uint** matrix, uint sizeR, uint sizeC, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < sizeR; i++) {
        cout << "[ ";
        for (uint j = 0; j < sizeC; j++) {
            cout << matrix[j][i] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}
void showMatrixT(int** matrix, uint sizeR, uint sizeC, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < sizeR; i++) {
        cout << "[ ";
        for (uint j = 0; j < sizeC; j++) {
            cout << matrix[j][i] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}
void showMatrixT(double** matrix, uint sizeR, uint sizeC, string name) {
    cout << name << " = " << endl;
    for (uint i = 0; i < sizeR; i++) {
        cout << "[ ";
        for (uint j = 0; j < sizeC; j++) {
            cout << matrix[j][i] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}

int** matMul(int** A, int** B, uint size) {
    int** C = new int*[size];
    for (uint i = 0; i < size; i++) {
        C[i] = new int[size];
    }
    for (uint i = 0; i < size; i++) {
        for (uint j = 0; j < size; j++) {
            for (uint k = 0; k < size; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}
double** matMul(double** A, double** B, uint size) {
    double** C = new double*[size];
    for (uint i = 0; i < size; i++) {
        C[i] = new double[size];
    }
    for (uint i = 0; i < size; i++) {
        for (uint j = 0; j < size; j++) {
            for (uint k = 0; k < size; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

int** outProd(int* vec1, int* vec2, uint size) {
    int** C = new int*[size];
    for (uint i = 0; i < size; i++) {
        C[i] = new int[size];
    }
    for (uint i = 0; i < size; i++) {
        for (uint j = 0; j < size; j++) {
            C[i][j] = vec1[i+1] * vec2[j];
        }
    }
    return C;
}
double** outProd(double* vec1, int* vec2, uint size) {
    double** C = new double*[size];
    for (uint i = 0; i < size; i++) {
        C[i] = new double[size];
    }
    for (uint i = 0; i < size; i++) {
        for (uint j = 0; j < size; j++) {
            C[i][j] = vec1[i+1] * vec2[j];
        }
    }
    return C;
}