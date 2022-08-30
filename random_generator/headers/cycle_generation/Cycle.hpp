#pragma once

#include "Matrix.hpp"
#include "update_generation/UpdateMatrix.hpp"

namespace randomgen {

  class Cycle {
  public:
    Matrix &getSlaterMatrix() { return matrix; }

    const Matrix &getSlaterMatrix() const { return matrix; }

    Matrix &getSlaterInverse() { return matrix_invt; }

    const Matrix &getSlaterInverseTransposed() const { return matrix_invt; }

    UpdateMatrix &getUpdateMatrix() { return update_array; }

    const UpdateMatrix &getUpdateMatrix() const { return update_array; }

  private:
    Matrix matrix;
    Matrix matrix_invt;
    UpdateMatrix update_array;
  };

}// namespace randomgen