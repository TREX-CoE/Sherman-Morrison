#pragma once

class Engine;

/**
 * @brief Interface for classes that execute a specific generation mode
 *
 * Subclassing this class allows one to easily add new custom generation modes to the engine.
 * (This class is purposefully simple to allow for easy subclassing, try to keep it this way)
 */
class EngineModeExecutor {
public:
  explicit EngineModeExecutor(Engine *parent_engine) : parent(parent_engine) {}

  virtual ~EngineModeExecutor() = default;

  /**
   * @brief Main method of the executor. This method is automatically called by the Engine class.
   * @return The return code of the execution, 0 on success, non-zero on error.
   */
  virtual int exec() = 0;

protected:
  /**
   * The engine that owns this executor
   * Can be used to fetch the program arguments, etc.
   */
  Engine *parent;
};
