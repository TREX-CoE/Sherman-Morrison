#pragma once
#include "UpdatePlan.hpp"
#include <iostream>
#include <vector>

namespace randomgen {

  class UpdatePlanner {
  public:
    virtual UpdatePlan make(uint32_t update_count, uint32_t split_count) const = 0;
  };


}// namespace randomgen
