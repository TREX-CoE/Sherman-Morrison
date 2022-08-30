
#pragma once

#include "EngineModeExecutor.hpp"

/**
 * @brief Generates a dataset with matrices of increasing size and increasing number of splits.
 * The total number of updates is fixed, and the number of splits is increasing.
 *
 * @TODO Add a method for customizing the executor's parameters
 */
class MatrixSizeExecutor : public EngineModeExecutor {
public:
  using EngineModeExecutor::EngineModeExecutor;

  int exec() override;
};