
#include "cycle_generation/SizedBasedCycleStream.hpp"

namespace randomgen {
  std::string SizedBasedCycleStream::getPathFor(const randomgen::Cycle &cycle) {
    // Save the cycle based on the size of its matrix + number of splits
    std::string res = "/";
    res += std::to_string(cycle.getSlaterMatrix().x) + "/";
    res += std::to_string(cycle.getUpdateMatrix().getSplitCount()) + "/";
    res += "cycle_" + std::to_string(last_unique_cycle_id);
    return res;
  }

}// namespace randomgen
