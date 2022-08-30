#pragma once

#include "UpdateDescriptor.hpp"
#include <iostream>
#include <vector>

namespace randomgen {

  /**
   * @brief Collection of UpdateDescriptors, describing all the updates in a cycle
   *
   * This class describes the updates contained in a cycle, as a set of UpdateDescriptors.
   * It can be fed to the UpdatePlanBuilder to generate a cycle.
   */
  class UpdatePlan {
  public:
    using Iterator = std::vector<UpdateDescriptor>::iterator;
    using ConstIterator = std::vector<UpdateDescriptor>::const_iterator;

    UpdatePlan() = default;

    UpdatePlan(uint32_t updates_count);

    Iterator begin() { return updates_descriptors.begin(); }

    Iterator end() { return updates_descriptors.end(); }

    ConstIterator begin() const { return updates_descriptors.begin(); }

    ConstIterator end() const { return updates_descriptors.end(); }

    UpdateDescriptor &operator[](uint32_t index) { return updates_descriptors[index]; }

    const UpdateDescriptor &operator[](uint32_t index) const { return updates_descriptors[index]; }

    UpdateDescriptor &at(uint32_t index) { return updates_descriptors.at(index); }

    const UpdateDescriptor &at(uint32_t index) const { return updates_descriptors.at(index); }

    uint32_t size() const { return updates_descriptors.size(); }

    uint32_t getSplitCount() const;

  private:
    std::vector<UpdateDescriptor> updates_descriptors;
  };

}// namespace randomgen
