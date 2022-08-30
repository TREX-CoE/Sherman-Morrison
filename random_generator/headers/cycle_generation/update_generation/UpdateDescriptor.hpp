#pragma once

#include <iostream>

namespace randomgen {


  /**
   * @brief Describes single update, including its target column, whether it splits or not...
   */
  class UpdateDescriptor {
  public:
    UpdateDescriptor() = default;

    /**
     * @brief Construct a new UpdateDescriptor object
     * @param splits If true, the update must splits (produce a singular matrix)
     * @param column_index The column to update
     * @param target_column_index The target column to use for the split (if the splits is False, this is ignored)
     *
     */
    UpdateDescriptor(bool splits, uint32_t column_index, uint32_t target_column_index)
        : splits(splits), column_index(column_index), target_column_index(target_column_index) {}

    bool isSplitting() const { return splits; }

    void setSplitting(bool new_split_value) { this->splits = new_split_value; }

    uint32_t getColumnIndex() const { return column_index; }

    void setColumnIndex(uint32_t new_column_index) { column_index = new_column_index; }

    uint32_t getTargetColumnIndex() const { return target_column_index; }

    void setTargetColumnIndex(uint32_t new_target_column_index) {
      this->target_column_index = new_target_column_index;
    }

  private:
    bool splits = false;
    uint32_t column_index = 0;
    uint32_t target_column_index = 0;
  };

}// namespace randomgen