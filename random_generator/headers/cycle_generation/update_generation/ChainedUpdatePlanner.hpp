#pragma once

#include <iostream>

#include "UpdatePlanner.hpp"

namespace randomgen {

  /**
   * @brief Generates an update plan using chained updates
   *
   * Updates are generated left-to-right, where each splitting update targets the column on its right.
   * By ensuring that the last update is non-splitting, we can guarantee that the matrix remains invertible.
   */
  class ChainedUpdatePlanBuilder : public UpdatePlanner {
  public:
    /**
     * @brief Produce a chain-update plan
     * @param update_count The total number of updates in the plan
     * @param split_count The number of splitting updates in the plan
     * @return A new plan with the given parameters
     */
    UpdatePlan make(uint32_t update_count, uint32_t split_count) const override;
  };

}// namespace randomgen
