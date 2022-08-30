#pragma once

#include "Cycle.hpp"
#include "Matrix.hpp"
#include "update_generation/UpdateGenerator.hpp"
#include "update_generation/UpdateMatrix.hpp"
#include <iostream>

namespace randomgen {

  class CycleGenerator {
  public:
    /**
     * @brief Cycle generator
     * @param x The size of the (squared) matrix
     * @param y_padding The column padding to add to the matrix
     * @param generator The generator used to generate the update vectors
     * Note that this generator must be thread-safe and re-entrant
     */
    CycleGenerator(size_t x, size_t y_padding,
                   std::shared_ptr<const UpdateGenerator> update_generator);

    void generateCycleMatrices(Cycle *res);

    Cycle make(size_t n_update, size_t n_splits);

  private:
    std::shared_ptr<const UpdateGenerator> update_generator;
    MatrixGenerator matrix_generator;
    size_t x, y_padding;

    Cycle attemptToGenerateCycle(size_t n_update, size_t n_splits);

    void finalizeCycle(Cycle &res) const;
  };
}// namespace randomgen
