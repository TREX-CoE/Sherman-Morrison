// Helpers.hpp
// Some usefull helper functions to support the Maponi algorithm.
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

template<typename T>
unsigned int getMaxIndex(T *vector, unsigned int size) {
    unsigned int i;
    unsigned int max = vector[0];
    unsigned int maxi = 0;
    for (i = 1; i < size; i++) {
        if (vector[i] > max) {
            max = vector[i];
            maxi = i;
        }
    }
    return maxi;
}

template<typename T>
void showScalar(T scalar, string name) {
    cout << name << " = " << scalar << endl << endl;
}

template<typename T>
void showVector(T *vector, unsigned int size, string name) {
    cout << name << " = " << endl;
    for (unsigned int i = 0; i < size; i++) {
        cout << "[ " << vector[i] << " ]" << endl;
    }
    cout << endl;
}

template<typename T>
void showMatrix(T *matrix, unsigned int M, string name) {
    cout << name << " = " << endl;
    for (unsigned int i = 0; i < M; i++) {
        cout << "[ ";
        for (unsigned int j = 0; j < M; j++) {
            cout << matrix[i*M+j] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}

template<typename T>
void showMatrixT(T **matrix, unsigned int size, string name) {
    cout << name << " = " << endl;
    for (unsigned int i = 0; i < size; i++) {
        cout << "[ ";
        for (unsigned int j = 0; j < size; j++) {
            cout << matrix[j][i] << " ";
        }
        cout << " ]" << endl;
    }
    cout << endl;
}

template<typename T>
T *matMul(T *A, T *B, unsigned int M) {
    T *C = new T[M*M];
    for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = 0; j < M; j++) {
            for (unsigned int k = 0; k < M; k++) {
                C[i*M+j] += A[i*M+k] * B[k*M+j];
            }
        }
    }
    return C;
}


template<typename T1, typename T2>
T1 *outProd(T1 *vec1, T2 *vec2, unsigned int M) {
    T1 *C = new T1[M*M];
    for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = 0; j < M; j++) {
            C[i*M+j] = vec1[i+1] * vec2[j];
        }
    }
    return C;
}

template<typename T>
T matDet(T **A, unsigned int M) {
    int det = 0, p, h, k, i, j;
    T **temp = new T*[M];
    for (int i = 0; i < M; i++) temp[i] = new T[M];
    if(M == 1) {
        return A[0][0];
    }
    else if(M == 2) {
        det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
        return det;
    }
    else {
        for(p = 0; p < M; p++) {
            h = 0;
            k = 0;
            for(i = 1; i < M; i++) {
                for( j = 0; j < M; j++) {
                    if(j == p) {
                        continue;
                    }
                    temp[h][k] = A[i][j];
                    k++;
                    if(k == M-1) {
                        h++;
                        k = 0;
                    }
                }
            }
            det = det + A[0][p] * pow(-1, p) * matDet(temp, M-1);
        }
        return det;
    }
    delete [] temp;
}
