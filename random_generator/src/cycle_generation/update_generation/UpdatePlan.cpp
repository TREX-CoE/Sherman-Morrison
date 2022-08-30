
#include "UpdatePlan.hpp"

namespace randomgen {

  UpdatePlan::UpdatePlan(uint32_t updates_count) { updates_descriptors.resize(updates_count); }

  uint32_t UpdatePlan::getSplitCount() const {
    uint32_t res = 0;
    for (const auto &update: updates_descriptors) { res += update.isSplitting(); }
    return res;
  }

}// namespace randomgen
