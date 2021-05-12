// SM_Helpers.hpp
// Some usefull helper functions to support the Maponi algorithm.
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

// #define DEBUG

#ifndef THRESHOLD
#define THRESHOLD 1e-3
#endif
double threshold();

void Switch(unsigned int *p, unsigned int l, unsigned int lbar);

void selectLargestDenominator(unsigned int l, unsigned int N_updates,
                              unsigned int *Updates_index, unsigned int *p,
                              double ***ylk);

template <typename T> void showScalar(T scalar, std::string name) {
  std::cout << name << " = " << scalar << std::endl << std::endl;
}

template <typename T>
void showVector(T *vector, unsigned int size, std::string name) {
  std::cout << name << " = " << std::endl;
  for (unsigned int i = 0; i < size; i++) {
    std::cout << "[ " << vector[i] << " ]" << std::endl;
  }
  std::cout << std::endl;
}

template <typename T>
void showMatrix(T *matrix, unsigned int M, std::string name) {
  std::cout.precision(17);
  std::cout << name << " = [" << std::endl;
  for (unsigned int i = 0; i < M; i++) {
    std::cout << "[";
    for (unsigned int j = 0; j < M; j++) {
      if (matrix[i * M + j] >= 0) {
        std::cout << "  " << matrix[i * M + j] << ",";
      } else {
        std::cout << " " << matrix[i * M + j] << ",";
      }
    }
    std::cout << " ]," << std::endl;
  }
  std::cout << "]" << std::endl;
  std::cout << std::endl;
}

template <typename T> T *transpose(T *A, unsigned int M) {
  T *B = new T[M * M];
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < M; j++) {
      B[i * M + j] = A[i + j * M];
    }
  }
  return B;
}

template <typename T> void matMul(T *A, T *B, T *C, unsigned int M) {
  memset(C, 0, M * M * sizeof(T));
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < M; j++) {
      for (unsigned int k = 0; k < M; k++) {
        C[i * M + j] += A[i * M + k] * B[k * M + j];
      }
    }
  }
}

template <typename T1, typename T2>
T1 *outProd(T1 *vec1, T2 *vec2, unsigned int M) {
  T1 *C = new T1[M * M];
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < M; j++) {
      C[i * M + j] = vec1[i + 1] * vec2[j];
    }
  }
  return C;
}

template <typename T> T matDet(T **A, unsigned int M) {
  int det = 0, p, h, k, i, j;
  T **temp = new T *[M];
  for (int i = 0; i < M; i++)
    temp[i] = new T[M];
  if (M == 1) {
    return A[0][0];
  } else if (M == 2) {
    det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
    return det;
  } else {
    for (p = 0; p < M; p++) {
      h = 0;
      k = 0;
      for (i = 1; i < M; i++) {
        for (j = 0; j < M; j++) {
          if (j == p) {
            continue;
          }
          temp[h][k] = A[i][j];
          k++;
          if (k == M - 1) {
            h++;
            k = 0;
          }
        }
      }
      det = det + A[0][p] * pow(-1, p) * matDet(temp, M - 1);
    }
    return det;
  }
  delete[] temp;
}

template <typename T> bool is_identity(T *A, unsigned int M, double tolerance) {
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < M; j++) {
      if (i == j && std::fabs(A[i * M + j] - 1) > tolerance)
        return false;
      if (i != j && std::fabs(A[i * M + j]) > tolerance)
        return false;
    }
  }
  return true;
}

template <typename T> T norm_max(T *A, unsigned int Dim) {
  T res = 0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T delta = A[i * Dim + j];
      delta = std::fabs(delta);
      if (delta > res) {
        res = delta;
      }
    }
  }
  return res;
}

template <typename T> T norm_frobenius2(T *A, unsigned int Dim) {
  T res = 0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T delta = A[i * Dim + j];
      res += delta * delta;
    }
  }
  return res;
}

template <typename T> T residual_max(T *A, unsigned int Dim) {
  T res = 0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T delta = A[i * Dim + j] - (i == j);
      delta = std::fabs(delta);
      if (delta > res) {
        res = delta;
      }
    }
  }
  return res;
}

template <typename T> T residual_frobenius2(T *A, unsigned int Dim) {
  T res = 0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T delta = A[i * Dim + j] - (i == j);
      res += delta * delta;
    }
  }
  return res;
}
