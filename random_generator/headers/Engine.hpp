#pragma once

#include <filesystem>
#include <iostream>

#include "EngineModeExecutor.hpp"

enum class EngineMode {
  kMatrixSize,
  kUpdateCount,
};

/**
 * @brief Main class responsible for starting the correct generation depending on the program arguments
 *
 * Adding a new generation mode is as simple as adding a new enum value to EngineMode and adapting the
 * exec() method
 */
class Engine {
public:
  Engine(int argc, char **argv);

  int exec();

  char *getExecutablePath() { return argv[0]; }

  const std::filesystem::path &getOutputPath() const { return dataset_output_path; }

  const char *const *getArgv() const { return argv; }

  int getArgc() const { return argc; }

private:
  char **argv;
  int argc;
  EngineMode mode;
  std::filesystem::path dataset_output_path;
};
