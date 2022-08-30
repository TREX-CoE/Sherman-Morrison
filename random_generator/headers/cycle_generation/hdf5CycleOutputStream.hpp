#pragma once


#include "Cycle.hpp"
#include <filesystem>
#include <hdf5.h>
#include <iostream>
#include <vector>


namespace randomgen {
  /**
   * @brief This stream allows one to output multiple cycles to a single HDF5 file.
   * The cycle storing algorithm is fixed, but this class can be subclassed to change where the cycles
   * are stored in the final datasets.
   *
   * By default, this class stores the cycles per number of updates and number of splits.
   */
  class hdf5CycleOutputStream {
  public:
    explicit hdf5CycleOutputStream(const std::filesystem::path &dataset_path);

    virtual ~hdf5CycleOutputStream();


    /**
     * @brief Output a cycle to the dataset.
     * @param ds The stream to output to.
     * @param cycle The cycle to output
     * @return The stream passed as parameter.
     */
    friend hdf5CycleOutputStream &operator<<(hdf5CycleOutputStream &ds, const Cycle &cycle);

  protected:
    /**
     * @brief Creates the given path inside the dataset (parents folders included)
     * @param path The path to create
     * @return the id of the newly created group
     * Freeing is the responsibility of the caller.
     */
    hid_t createGroup(const std::string &path);

    hid_t file_id = -1;

    // Each cycle is represented by a unique id
    // This member is used to keep track of the current id
    size_t last_unique_cycle_id = 0;

  private:
    /**
     * Get the storage path for a given cycle.
     *
     * By default, this is /n_updates/n_splits
     * @param cycle The cycle that will be stored
     * @return The path at which the cycle must be stored.
     */
    virtual std::string getPathFor(const Cycle &cycle);
  };
}// namespace randomgen
