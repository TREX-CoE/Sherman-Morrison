
#include "cycle_generation/Matrix.hpp"
#include <lapacke.h>

namespace randomgen {

  Matrix::Matrix(size_t x, size_t y) : x(x), y(y) { array = std::make_unique<double[]>(x * y); }

  Matrix::Matrix(const Matrix &other) : Matrix(other.x, other.y) {
    std::copy(other.array.get(), other.array.get() + x * y, array.get());
    determinant = other.determinant;
  }

  Matrix &Matrix::operator=(const Matrix &other) {
    if (this == &other) { return *this; }
    // Realloc a matrix if the size is different

    if (x * y != other.x * other.y) { array = std::make_unique<double[]>(other.x * other.y); }

    x = other.x;
    y = other.y;

    std::copy(other.array.get(), other.array.get() + x * y, array.get());
    return *this;
  }

  Matrix Matrix::inverse() {
    // Return the inverse of this matrix using lapacke dgetrf + dgetri
    Matrix inv = *this;
    std::vector<int> ipiv(x);
    auto info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, x, y, inv.array.get(), x, ipiv.data());
    if (info != 0) { throw std::runtime_error("Matrix inversion failed"); }
    auto info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, x, inv.array.get(), x, ipiv.data());
    if (info2 != 0) { throw std::runtime_error("Matrix inversion failed"); }

    inv.computeDeterminant();
    return inv;
  }

  Matrix Matrix::transpose() {
    Matrix res(x, y);

    for (size_t i = 0; i < x; i++) {
      for (size_t j = 0; j < y; j++) { res.array[j * x + i] = array[i * y + j]; }
    }
    res.computeDeterminant();

    return res;
  }

  std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    for (size_t i = 0; i < m.x; i++) {
      for (size_t j = 0; j < m.y; j++) { os << m.array[i * m.y + j] << " "; }
      os << std::endl;
    }
    return os;
  }

  void Matrix::pad(size_t x_pad, size_t y_pad) {
    if (x_pad == 0 && y_pad == 0) { return; }

    size_t new_size = (x + x_pad) * (y + y_pad);
    std::unique_ptr<double[]> padded_data = std::make_unique<double[]>(new_size);

    // Set the padding to 0
    std::fill(padded_data.get(), padded_data.get() + (new_size), 0);

    // Copy the old data to the new data
    for (int i = 0; i < x; i++) {
      for (int j = 0; j < y; j++) { padded_data[i * (y + y_pad) + j] = array[i * y + j]; }
    }
    x += x_pad;
    y += y_pad;
    this->array = std::move(padded_data);
  }

  double Matrix::computeDeterminant() {
    // Returns the determinant of the input matrix using LAPACK
    // Copy the matrix as lapack modifies it
    std::unique_ptr<double[]> copy = std::make_unique<double[]>(x * y);
    std::copy(array.get(), array.get() + (x * y), copy.get());

    std::vector<int> pivot(x);
    auto info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, x, y, copy.get(), x, pivot.data());
    if (info != 0) {
      std::cout << "HERE!" << std::endl;
      throw std::runtime_error("Matrix inversion failed");
    }

    double det = 1;
    for (int i = 0; i < y; ++i) { det *= array[i * x + i]; }
    determinant = det;
    return det;
  }

  double Matrix::getLastKnownDeterminant() const { return determinant; }


  Matrix MatrixGenerator::operator()(size_t matrix_size, double vmin, double vmax) {
    Matrix res{matrix_size, matrix_size};
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(vmin, vmax);
    std::generate(res.array.get(), res.array.get() + (matrix_size * matrix_size),
                  [&]() { return dis(gen); });
    res.computeDeterminant();
    return res;
  }

}// namespace randomgen