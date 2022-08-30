
#include "ChainedUpdatePlanner.hpp"
#include <algorithm>
#include <random>

namespace randomgen {


  namespace {

    /**
     * @brief Iterate on a plan and reorders the split to ensure the last update does not split
     * @param plan
     */
    void reorderSplits(UpdatePlan *plan) {

      // Find the first non-splitting updates
      size_t first_non_split = -1;
      for (size_t i = 0; i < plan->size(); ++i) {
        if (not plan->operator[](i).isSplitting()) {
          first_non_split = i;
          break;
        }
      }

      if (first_non_split == -1) {
        throw std::runtime_error("Reordering Failure: No non-splitting update available for "
                                 "reordering was found in the plan");
      }

      // Swap the first non-splitting update with the last splitting update
      plan->operator[](plan->size() - 1).setSplitting(false);
      plan->operator[](first_non_split).setSplitting(true);
    }

    /**
     * @brief Check the order of the splitting updates in the plan, and maybe reorder them if the last update is
     * a splitting update
     *
     * @param plan
     * @param splits
     */
    void maybeReorderSplits(UpdatePlan *plan, const size_t splits) {
      if (plan->size() < 2 or splits == 0) { return; }

      // Ensure that the cycle is broken by the last update
      // Should look something like [..., splits, doesn't splits]
      bool invalid_chain = plan->at(plan->size() - 1).isSplitting();
      if (not invalid_chain) { return; }

      // If the chain is invalid, reorder the splits
      reorderSplits(plan);
    }

    /**
     * @brief Make a plan where updates are generated left-to-right, and splitting updates are based on the right column
     * This allows the colinearity of the updates to be broken by the next update, effectively ensuring that the final matrix
     * is not singular
     *
     * @param plan
     * @param splits
     */
    void makeChainedUpdatePlan(UpdatePlan *plan, const size_t n_update, const size_t splits) {

      std::vector<bool> split_order(plan->size(), false);
      for (size_t i = 0; i < splits; ++i) { split_order[i] = true; }

      std::shuffle(split_order.begin(), split_order.end(), std::mt19937{std::random_device{}()});

      for (size_t i = 0; i < plan->size(); ++i) {
        auto &update_desc = plan->at(i);
        update_desc.setSplitting(split_order[i]);
        // Always split by
        update_desc.setTargetColumnIndex(i + 1);
        update_desc.setColumnIndex(i);
      }

      maybeReorderSplits(plan, splits);
    }

  }// namespace

  UpdatePlan ChainedUpdatePlanBuilder::make(uint32_t update_count, uint32_t split_count) const {
    if (split_count >= update_count)
      throw std::runtime_error("Cannot have more splits than updates");

    UpdatePlan res(update_count);
    makeChainedUpdatePlan(&res, update_count, split_count);
    return res;
  }

}// namespace randomgen
