// Helpers.hpp
// Some usefull helper functions to support the Maponi algorithm.
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#ifdef MKL
#include <mkl_lapacke.h>
#endif

// #define DEBUG
#ifndef THRESHOLD
#define THRESHOLD 1e-3
#endif

double threshold();

void Switch(unsigned int *p, unsigned int l, unsigned int lbar);

void selectLargestDenominator(unsigned int l, unsigned int N_updates,
                              unsigned int *Updates_index, unsigned int *p,
                              double ***ylk);

#ifdef MKL
lapack_int inverse(double *A, unsigned n);
#endif

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

template <typename T>
void showMatrix2(T *matrix, unsigned int M, unsigned int N, std::string name) {
  std::cout.precision(17);
  std::cout << name << " = [" << std::endl;
  for (unsigned int i = 0; i < M; i++) {
    std::cout << "[";
    for (unsigned int j = 0; j < N; j++) {
      if (matrix[i * N + j] >= 0) {
        std::cout << "  " << matrix[i * N + j] << ",";
      } else {
        std::cout << " " << matrix[i * N + j] << ",";
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

template <typename T1, typename T2, typename T3>
void matMul2(T1 *A, T2 *B, T3 *C, unsigned int M, unsigned int N, unsigned int P) {
  for(unsigned int i = 0; i < M; i++) {
    for(unsigned int j = 0; j < P; j++) {
      C[i * P + j] = 0;
      for(unsigned int k = 0; k < N; k++) {
        C[i * P + j] += A[i * N + k] * B[k * P + j];
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

// // This flat version doesn't work. Get's stuck in an infinite recursion loop.
// template <typename T> T determinant(T *A, unsigned int M) {
//   std::cout << "determinant() called..." << std::endl;
//   T det = 0;
//   int p, h, k, i, j;
//   T *temp = new T[M * M];
//   if (M == 1) {
//     return A[0];
//   } else if (M == 2) {
//     det = (A[0] * A[3] - A[1] * A[2]);
//     return det;
//   } else {
//     for (p = 0; p < M; p++) {
//       h = 0;
//       k = 0;
//       for (i = 1; i < M; i++) {
//         for (j = 0; j < M; j++) {
//           if (j == p) {
//             continue;
//           }
//           temp[h * M + k] = A[i * M + j];
//           k++;
//           if (k == M - 1) {
//             h++;
//             k = 0;
//           }
//         }
//       }
//       det = det + A[p] * pow(-1, p) * determinant(temp, M - 1);
//     }
//     return det;
//   }
//   delete temp;
// }

// // This version also gets stuck in a recursion loop
// template <typename T> T determinant(T **A, unsigned int M) {
//   int p, h, k, i, j;
//   T det = 0;
//   T **temp = new T *[M];
//   for (int i = 0; i < M; i++) {
//     temp[i] = new T[M];
//   }
//   if (M == 1) {
//     return A[0][0];
//   } else if (M == 2) {
//     det = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
//     return det;
//   } else {
//     for (p = 0; p < M; p++) {
//       h = 0;
//       k = 0;
//       for (i = 1; i < M; i++) {
//         for (j = 0; j < M; j++) {
//           if (j == p) {
//             continue;
//           }
//           temp[h][k] = A[i][j];
//           k++;
//           if (k == M - 1) {
//             h++;
//             k = 0;
//           }
//         }
//       }
//       det = det + A[0][p] * pow(-1, p) * determinant(temp, M - 1);
//     }
//     return det;
//   }
//   delete[] temp;
// }

template <typename T> bool is_identity(T *A, unsigned int M, double tolerance) {
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < M; j++) {
      if (i == j && std::fabs(A[i * M + j] - 1) > tolerance) {
        return false;
      }
      if (i != j && std::fabs(A[i * M + j]) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

template <typename T>
bool is_identity2(T *A, unsigned int M, double tolerance) {
  double det = determinant(A, M);
  if (det - 1 > tolerance) {
    return false;
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


template <typename T> T residual2(T * A, unsigned int Dim) {
  double res = 0.0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T delta = (A[i * Dim + j] - (i == j));
      res += delta*delta;
    }
  }
  return res;
}

// Computes the condition number of A using a previously computed B=A^{-1}
template <typename T> T condition1(T *A, T *B, unsigned int Dim) {
  T resA = 0;
  T resB = 0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T deltaA = A[i * Dim + j];
      T deltaB = B[i * Dim + j];
      resA += deltaA * deltaA;
      resB += deltaB * deltaB;
    }
  }
  return sqrt(resA * resB);
}

#ifdef MKL
// Computes the condition number of A by first inverting it with LAPACKE
template <typename T> T condition2(T *A, unsigned int Dim) {
  T B[Dim * Dim];
  std::memcpy(B, A, Dim * Dim * sizeof(T));
  inverse(B, Dim);
  T resA = 0;
  T resB = 0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      T deltaA = A[i * Dim + j];
      T deltaB = B[i * Dim + j];
      resA += deltaA * deltaA;
      resB += deltaB * deltaB;
    }
  }
  return sqrt(resA) * sqrt(resB);
}
#endif
