#pragma once
#include "hdf5CycleOutputStream.hpp"


namespace randomgen {

  /**
   * @brief A stream that stores cycles based on the matrices size
   */
  class SizedBasedCycleStream : public hdf5CycleOutputStream {
  public:
    using hdf5CycleOutputStream::hdf5CycleOutputStream;

  private:
    /**
     * @brief Returns the path as /matrix_size/splits_count/cycle_id
     *
     * @param cycle
     * @return
     */
    std::string getPathFor(const Cycle &cycle) override;
  };

}// namespace randomgen
