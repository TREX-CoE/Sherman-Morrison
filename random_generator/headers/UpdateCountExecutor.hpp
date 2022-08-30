#pragma once

#include "EngineModeExecutor.hpp"

/**
 * @brief Generates a dataset where cycles contain a varying number of updates.
 * The number of splitting updates varies between 2 and the number of updates.
 *
 * @TODO Add a method for customizing the executor's parameters
 */
class UpdateCountExecutor : public EngineModeExecutor {
public:
  using EngineModeExecutor::EngineModeExecutor;

  int exec() override;
};