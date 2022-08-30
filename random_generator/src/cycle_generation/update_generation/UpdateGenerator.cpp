
#include "UpdateGenerator.hpp"

#include <utility>

namespace randomgen {

  UpdateGenerator::UpdateGenerator(std::shared_ptr<const UpdatePlanner> update_planner,
                                   std::shared_ptr<const UpdateMatrixGenerator> update_generator)
      : planner(std::move(update_planner)), generator(std::move(update_generator)) {}


  UpdatePlan UpdateGenerator::plan(uint32_t update_count, uint32_t split_count) const {
    return planner->make(update_count, split_count);
  }

  UpdateMatrix UpdateGenerator::make(const Matrix &slater_matrix, uint32_t y_padding,
                                     const UpdatePlan &update_plan) const {
    return generator->make(slater_matrix, y_padding, update_plan);
  }


  UpdateMatrix UpdateGenerator::make(const Matrix &slater_matrix_t, uint32_t y_padding,
                                     uint32_t update_count, uint32_t split_count) const {
    auto new_plan = plan(update_count, split_count);
    return make(slater_matrix_t, y_padding, new_plan);
  }


}// namespace randomgen