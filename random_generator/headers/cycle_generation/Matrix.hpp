#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

namespace randomgen {
  struct Matrix {
    Matrix() = default;
    Matrix(size_t x, size_t y);

    Matrix(const Matrix &other);
    Matrix &operator=(const Matrix &other);

    Matrix(Matrix &&other) = default;
    Matrix &operator=(Matrix &&other) = default;

    Matrix inverse();
    Matrix transpose();

    friend std::ostream &operator<<(std::ostream &os, const Matrix &m);

    double *data() { return array.get(); }

    const double *data() const { return array.get(); }

    void pad(size_t x_pad, size_t y_pad);
    double computeDeterminant();
    double getLastKnownDeterminant() const;

    uint32_t x = 0, y = 0;
    double determinant;
    std::unique_ptr<double[]> array;
  };

  class MatrixGenerator {
  public:
    Matrix operator()(size_t matrix_size, double vmin = -3, double vmax = 3);
  };
}// namespace randomgen