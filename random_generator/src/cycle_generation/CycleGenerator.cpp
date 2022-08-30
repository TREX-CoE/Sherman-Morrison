

#include "cycle_generation/CycleGenerator.hpp"

#include <utility>

namespace randomgen {

  CycleGenerator::CycleGenerator(size_t x, size_t y_padding,
                                 std::shared_ptr<const UpdateGenerator> update_generator)
      : update_generator(std::move(update_generator)), x(x), y_padding(y_padding) {}


  void CycleGenerator::generateCycleMatrices(Cycle *res) {
    // Generate random matrices until we find one with a given determinant
    do {
      res->getSlaterMatrix() = matrix_generator(x);
    } while (res->getSlaterMatrix().getLastKnownDeterminant() < 1e-5);

    res->getSlaterInverse() = res->getSlaterMatrix().inverse();
  }


  Cycle CycleGenerator::attemptToGenerateCycle(size_t n_update, size_t n_splits) {
    Cycle res{};
    size_t attempt = 0;
    constexpr size_t kMaxAttempts = 20;

    auto plan = update_generator->plan(n_update, n_splits);

    while (true) {
      attempt++;
      generateCycleMatrices(&res);

      // Try to find a valid update vector
      // If we can't, start from scratch
      try {
        res.getUpdateMatrix() =
                update_generator->make(res.getSlaterMatrix().transpose(), y_padding, plan);
        break;
      } catch (UpdateMatrixGenerator::FailedUpdateGeneration &e) {
        // Abort after a certain number of attempts
        if (attempt > kMaxAttempts) {
          throw std::runtime_error("Could not generate a valid cycle in " +
                                   std::to_string(attempt) + " attempts");
        }
        continue;
      }
    }
    return res;
  }

  void CycleGenerator::finalizeCycle(Cycle &res) const {
    // The inverse in the dataset must be transposed
    // this allows the code of Scherman-morrison to be vectorized since column-updates becomes row-updates
    res.getSlaterInverse() = res.getSlaterInverse().transpose();
    res.getSlaterMatrix().pad(0, y_padding);
    res.getSlaterInverse().pad(0, y_padding);

    for (uint32_t i = 0; i < res.getUpdateMatrix().getUpdateCount(); i++) {
      res.getUpdateMatrix().getUpdateIndices()[i] += 1;
    }
  }

  Cycle CycleGenerator::make(size_t n_update, size_t n_splits) {

    Cycle res = attemptToGenerateCycle(n_update, n_splits);
    finalizeCycle(res);

    return res;
  }


}// namespace randomgen