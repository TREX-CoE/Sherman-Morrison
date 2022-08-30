#pragma once

#include "UpdateMatrixGenerator.hpp"
#include "UpdatePlanner.hpp"
#include <iostream>

namespace randomgen {

  /**
   * @brief Main class for the generator module, capable of planning and generating updates
   *
   * This class is a utility class that combines a planner and a generator, into a single entry point for the module.
   * It provides simplicity for the end user, while advanced user may want to use the planner and generator separately.
   */
  class UpdateGenerator {
  public:
    /**
     * @brief Builds a new generator from a planner and an update generator
     *
     * Takes unique_ptr to ensure thread-safety and re-entry
     * @param update_planner The planner to use for generating plans
     * @param update_generator The generator used to build an UpdateMatrix from a plan
     */
    UpdateGenerator(std::shared_ptr<const UpdatePlanner> update_planner,
                    std::shared_ptr<const UpdateMatrixGenerator> update_generator);

    /**
     * @brief Generate a new update plan to be fed to the generator at a later time
     * Allows one to generate a plan that can be reused multiple times.
     *
     * This function is thread-safe and re-entrant.
     * @param update_count The number of updates to plan for
     * @param split_count The number of splitting updates inside the plan
     * @return A new plan with the given parameters
     */
    UpdatePlan plan(uint32_t update_count, uint32_t split_count) const;

    /**
     * @brief Run the generator on the given update plan
     *
     * This function is thread-safe and re-entrant.
     * @param slater_matrix The slater matrix to use for the update generation
     * @param y_padding The padding to add to the update matrix.
     * Initially, the updates have the same length as a column of the slater matrix
     *
     * @param update_plan
     * @return
     */
    UpdateMatrix make(const Matrix &slater_matrix, uint32_t y_padding,
                      const UpdatePlan &update_plan) const;

    /**
     * @brief Plans and generate an update matrix from scratch.
     *
     * Note that a new plan is created each time this function is called. For maximum performance,
     * consider using the plan() function to generate a plan in situations where the same plan can be reused
     * multiple times.
     *
     * This function is thread-safe and re-entrant.
     * @param slater_matrix_t
     * @param y_padding
     * @param update_count
     * @param split_count
     * @return
     */
    UpdateMatrix make(const Matrix &slater_matrix_t, uint32_t y_padding, uint32_t update_count,
                      uint32_t split_count) const;

  private:
    std::shared_ptr<const UpdatePlanner> planner;
    std::shared_ptr<const UpdateMatrixGenerator> generator;
  };

}// namespace randomgen