// SM-MaponiA3.hpp
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

uint getMaxIndex(double *arr, uint size);
template<typename T>void showScalar(T scalar, string name);
template<typename T>void showVector(T *vector, uint size, string name);
template<typename T>void showMatrix(T **matrix, uint size, string name);
template<typename T>void showMatrixT(T **matrix, uint size, string name);
template<typename T>T **matMul(T **A, T **B, uint size);
template<typename T1, typename T2>T1 **outProd(T1 *vec1, T2 *vec2, uint size);
template<typename T>T matDet(T **A, int M);
void Sherman_Morrison(int **Slater0, double **Slater_inv, uint *Dim, uint *N_updates, int **Updates, uint *Updates_index);
