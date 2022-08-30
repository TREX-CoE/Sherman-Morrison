
#include "UpdateMatrixGenerator.hpp"

namespace randomgen {

  namespace {


    std::vector<double> tryGenerateNonSplittingUpdate(const UpdateDescriptor &descriptor,
                                                      const Matrix &current_matrix) {

      size_t kMaxAttempts = 100;
      size_t attempts = 0;

      Matrix copy = current_matrix;
      std::vector<double> update(current_matrix.x);

      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(0, 5);

      while (true) {
        // Generate a new update vector
        std::generate(update.begin(), update.end(), [&]() { return dis(gen); });

        // Check if the update vector is a valid update
        for (int i = 0; i < current_matrix.x; i++) {
          copy.array[i * current_matrix.y + descriptor.getColumnIndex()] += update[i];
        }

        // After the update is applied, the determinant of the matrix should be > 10e-3
        double absolute_det = std::abs(copy.computeDeterminant());
        if (absolute_det > 1e-3) { break; }

        // Reverse the update
        for (int i = 0; i < current_matrix.x; i++) {
          copy.array[i * current_matrix.y + descriptor.getColumnIndex()] -= update[i];
        }

        attempts++;
        if (attempts > kMaxAttempts) {
          throw UpdateMatrixGenerator::FailedUpdateGeneration(attempts);
        }
      }

      return update;
    }

    std::vector<double> &addRandomNoise(std::vector<double> &update, const double noise_magnitude) {
      // We add a noise of magnitude 1e-6 to the update vector so the resulting columns are not perfectly collinear
      // Which ensures the matrix stays invertible.

      // First generate a random vector
      std::vector<double> noise(update.size());
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(0, 100);

      std::generate(noise.begin(), noise.end(), [&]() { return dis(gen); });

      // Compute the norm
      double norm = 0;
      for (size_t i = 0; i < update.size(); i++) { norm += noise[i] * noise[i]; }

      // Coefficient to have a magnitude of 1e-6
      // No need to update the noise vector, just multiply it in place
      norm = noise_magnitude / sqrt(norm);
      for (size_t i = 0; i < update.size(); i++) { update[i] += noise[i] * norm; }
      return update;
    }

    std::vector<double> generate_splitting_update(const UpdateDescriptor &descriptor,
                                                  const Matrix &m, const double noise_magnitude) {

      size_t target_index = descriptor.getTargetColumnIndex();

      constexpr double lambda = 1;
      std::vector<double> update(m.x);
      for (size_t i = 0; i < m.x; i++) {
        update[i] = (lambda * m.array[i * m.y + target_index])       // Target column
                    - m.array[i * m.y + descriptor.getColumnIndex()];// Current column
      }

      // Add a small random noise to the update vector to break the singularity
      update = addRandomNoise(update, noise_magnitude);

      return update;
    }

    /**
     * The updates may be padded with extra zeros, so we need to check the generated update size
     */
    void maybePadUpdate(std::vector<double> &update_vec, uint32_t padded_size) {
      if (padded_size == update_vec.size()) return;

      if (padded_size < update_vec.size())
        throw std::runtime_error("Padded update size is smaller than the original update size");

      size_t old_size = update_vec.size();
      update_vec.resize(padded_size);
      for (size_t i = old_size; i < padded_size; i++) { update_vec[i] = 0; }
    }

    void updateSlaterMatrix(
            Matrix &slater_matrix, const UpdateDescriptor &descriptor,
            const std::vector<double> &
                    update_vec) {// We found the correct update vector, so return it + update the current matrix
      for (int i = 0; i < slater_matrix.y; i++) {
        slater_matrix.array[i * slater_matrix.y + descriptor.getColumnIndex()] += update_vec[i];
      }
    }

    std::vector<double> makePotentialUpdate(const UpdateDescriptor &descriptor,
                                            const Matrix &slater_matrix, uint32_t padded_size,
                                            double noise_magnitude) {

      std::vector<double> res;
      if (descriptor.isSplitting()) {
        res = generate_splitting_update(descriptor, slater_matrix, noise_magnitude);
      } else {
        res = tryGenerateNonSplittingUpdate(descriptor, slater_matrix);
      }

      maybePadUpdate(res, padded_size);
      return res;
    }

  }// namespace

  // Alias for the exception throw when the generation fails to find a correct update.
  // This exception is not make for runtime error !
  using FailedUpdateGeneration = UpdateMatrixGenerator::FailedUpdateGeneration;

  FailedUpdateGeneration::FailedUpdateGeneration(size_t attempts)
      : std::runtime_error("Failed to generate a valid update after " + std::to_string(attempts) +
                           " attempts"),
        attempted_generation(attempts) {}

  Matrix UpdateMatrixGenerator::iterateOnPlan(const Matrix &m, UpdateMatrix *update_matrix,
                                              const UpdatePlan &plan) const {
    // This matrix represents the matrix after the nth update
    // At the beginning, it is the same as the input matrix
    Matrix current_slater_matrix = m;


    // Iterate on the plan, generating fitting updates
    for (uint32_t rank = 0; const auto &descriptor: plan) {

      // Randomly generate a potential update
      const auto padded_size = update_matrix->getUpdateLength();
      auto update_vec =
              makePotentialUpdate(descriptor, current_slater_matrix, padded_size, noise_magnitude);
      update_matrix->setUpdate(rank, update_vec, descriptor.getColumnIndex());
      updateSlaterMatrix(current_slater_matrix, descriptor, update_vec);

      rank++;
    }
    return current_slater_matrix;
  }

  bool UpdateMatrixGenerator::tryGenerateUpdates(UpdateMatrix *update_matrix,
                                                 const UpdatePlan &plan, const Matrix &m) const {

    Matrix resulting_matrix = iterateOnPlan(m, update_matrix, plan);

    // Check if the final matrix is invertible or not
    const double det = resulting_matrix.getLastKnownDeterminant();
    if (fabs(det) < determinant_threshold) { return false; }
    return true;
  }

  void UpdateMatrixGenerator::attemptToGenerateUpdates(const UpdatePlan &plan, UpdateMatrix *res,
                                                       const Matrix &m) const {
    size_t attempts = 0;
    constexpr size_t max_attempts = 10;

    // Attempt to generate all the updates until we find a valid series, or we reached the maximum number of attempts
    while (true) {
      bool success = tryGenerateUpdates(res, plan, m);

      if (success) { return; }

      attempts++;
      if (attempts >= max_attempts) { throw FailedUpdateGeneration(attempts); }
    }
  }

  UpdateMatrix UpdateMatrixGenerator::make(const Matrix &m, size_t y_padding,
                                           const UpdatePlan &plan) const {

    size_t update_length = m.x + y_padding;

    UpdateMatrix res(plan.size(), plan.getSplitCount(), update_length);
    attemptToGenerateUpdates(plan, &res, m);

    return res;
  }

}// namespace randomgen
