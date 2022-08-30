#pragma once

#include "UpdateMatrix.hpp"
#include "UpdatePlan.hpp"
#include "cycle_generation/Matrix.hpp"
#include <iostream>
#include <span>

namespace randomgen {

  /**
   * @brief Cycle updates generator that handles splitting and normal cycles
   */
  class UpdateMatrixGenerator {
  public:
    /**
     * @brief Exception thrown when the update generator fails to generate an update after a set number of attempts
     */
    class FailedUpdateGeneration : std::runtime_error {
    public:
      explicit FailedUpdateGeneration(size_t attempts);

      /**
       * @brief Get the number of attempts that were made to generate an update before failing
       * @return
       */
      size_t getAttemptedGenerationCount() const { return attempted_generation; }

    private:
      const size_t attempted_generation;
    };

    /**
     * @brief Construct a new generator with the given thresholds
     * @param noise_magnitude The magnitude of the noise added to splitting updates.
     * The noise is added so that the final matrix is not singular.
     * If the noise is too small, Scherman-Morrison will fail, if it is too large, the update may not trigger any splits.
     * Tested with a value of 1e-6
     * @param determinant_threshold
     */
    explicit UpdateMatrixGenerator(double noise_magnitude, double determinant_threshold = 1e-6)
        : noise_magnitude(noise_magnitude), determinant_threshold(determinant_threshold) {}


    /**
     * @brief Main function of the generator, takes a matrix as an input, and generate a series of update vectors
     * Note that this method is thread-safe and re-entrant.
     * @param slater_matrix The reference matrix
     * @param update_plan The number of updates to generate
     * @return An array of update vectors
     */
    UpdateMatrix make(const Matrix &slater_matrix, size_t y_padding,
                      const UpdatePlan &update_plan) const;

  private:
    const double noise_magnitude;
    const double determinant_threshold;

    Matrix iterateOnPlan(const Matrix &m, UpdateMatrix *update_matrix,
                         const UpdatePlan &plan) const;

    void attemptToGenerateUpdates(const UpdatePlan &plan, UpdateMatrix *res, const Matrix &m) const;

    bool tryGenerateUpdates(UpdateMatrix *update_matrix, const UpdatePlan &plan,
                            const Matrix &m) const;
  };


}// namespace randomgen