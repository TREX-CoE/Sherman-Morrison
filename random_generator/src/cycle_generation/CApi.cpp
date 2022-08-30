#include "CApi.h"
#include "ChainedUpdatePlanner.hpp"
#include "CycleGenerator.hpp"
#include <iostream>
#include <memory>

namespace {

  void convertMatrix(randomgen::Matrix source, double **destination) {
    *destination = (double *) malloc(sizeof(double) * source.x * source.y);
    std::copy(source.data(), source.data() + source.x * source.y, *destination);
  }

  void convertMatrices(randomgen::Cycle &cpp_cycle, Cycle *res) {// Copy the slater matrix
    // Copy the inverse slater matrix
    convertMatrix(cpp_cycle.getSlaterMatrix(), &res->slater_matrix);
    convertMatrix(cpp_cycle.getSlaterInverseTransposed(), &res->slater_inverse_t);
  }


  void convertUpdates(randomgen::Cycle &cpp_cycle, Cycle *res) {
    auto &update_matrix = cpp_cycle.getUpdateMatrix();
    auto raw_updates = update_matrix.getRawUpdates();
    res->updates = (double *) malloc(sizeof(double) * update_matrix.getUpdateCount() *
                                     update_matrix.getUpdateLength());

    std::copy(raw_updates,
              raw_updates + update_matrix.getUpdateCount() * update_matrix.getUpdateLength(),
              res->updates);

    res->col_update_index = (uint32_t *) malloc(sizeof(uint32_t) * update_matrix.getUpdateCount());
    std::copy(update_matrix.getUpdateIndices(),
              update_matrix.getUpdateIndices() + update_matrix.getUpdateCount(),
              res->col_update_index);
  }


  Cycle *convertCPPCycleToC(randomgen::Cycle &cpp_cycle) {
    auto *res = (Cycle *) malloc(sizeof(Cycle));
    res->lds = cpp_cycle.getSlaterMatrix().y;
    res->dim = cpp_cycle.getSlaterMatrix().x;
    res->n_update = cpp_cycle.getUpdateMatrix().getUpdateCount();
    res->determinant = cpp_cycle.getSlaterMatrix().getLastKnownDeterminant();

    convertMatrices(cpp_cycle, res);
    convertUpdates(cpp_cycle, res);

    return res;
  }

}// namespace

// Those functions must have C linkage
extern "C" {

Cycle *generateRandomCycle(uint32_t dim, uint32_t lds, uint32_t n_update, uint32_t n_splits,
                           double determinant_threshold, double splitting_update_noise_magnitude) {

  if (lds < dim) {
    std::cerr << "generateRandomCycle: lds < dim" << std::endl;
    return nullptr;
  }

  auto y_padding = lds - dim;

  // We need to build a generator at every call, which is not optimal
  auto planner = std::make_unique<randomgen::ChainedUpdatePlanBuilder>();
  auto matrix_generator = std::make_unique<randomgen::UpdateMatrixGenerator>(
          splitting_update_noise_magnitude, determinant_threshold);
  auto update_generator = std::make_shared<randomgen::UpdateGenerator>(std::move(planner),
                                                                       std::move(matrix_generator));


  // This will generate a random cycle
  // However, those class use smart pointers (that are private to the library)
  // So we can't return them directly.
  // This requires the additional step of converting the C++ class to a C struct.
  // Which has a cost memory/performance wise
  randomgen::CycleGenerator generator(dim, y_padding, update_generator);
  randomgen::Cycle cpp_cycle = generator.make(n_update, n_splits);
  auto converted = convertCPPCycleToC(cpp_cycle);
  return converted;
}

void freeCycle(Cycle **cycle) {
  if (cycle == nullptr or *cycle == nullptr) {
    std::cerr << "freeCycle: cycle is nullptr !" << std::endl;
    return;
  }


  auto c = *cycle;
  free(c->slater_matrix);
  free(c->slater_inverse_t);
  free(c->updates);
  free(c->col_update_index);
  free(c);
  *cycle = nullptr;
}
}
