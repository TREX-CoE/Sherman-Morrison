
#include "UpdateMatrix.hpp"


namespace randomgen {

  UpdateMatrix::UpdateMatrix(const UpdateMatrix &other)
      : UpdateMatrix(other.n_update, other.n_splits, other.update_length) {}

  UpdateMatrix::UpdateMatrix(uint32_t n_update, uint32_t n_splits, uint32_t update_length)
      : n_update(n_update), n_splits(n_splits), update_length(update_length) {
    matrix = std::make_unique<double[]>(n_update * update_length);
    update_index = std::make_unique<uint32_t[]>(n_update);
  }

  UpdateMatrix &UpdateMatrix::operator=(const UpdateMatrix &other) {
    if (this == &other) { return *this; }

    // Reallocate arrays if needed
    if (n_update * update_length != other.n_update * other.update_length) {
      matrix = std::make_unique<double[]>(other.n_update * other.update_length);
    }

    if (n_update != other.n_update) { update_index = std::make_unique<uint32_t[]>(other.n_update); }

    n_update = other.n_update;
    n_splits = other.n_splits;
    update_length = other.update_length;

    std::copy(other.matrix.get(), other.matrix.get() + n_update * update_length, matrix.get());
    std::copy(other.update_index.get(), other.update_index.get() + n_update, update_index.get());
    return *this;
  }

  void UpdateMatrix::setUpdate(uint32_t update_rank, const std::span<const double> &update,
                               uint32_t column_index) {
    // Fetch the view where the update must be stored inside the update matrix
    auto update_span = getUpdateSpan(update_rank);

    if (update.size() > update_span.size()) {
      throw std::runtime_error("Update vector is too large");
    }

    std::copy(update.begin(), update.end(), update_span.begin());
    update_index.get()[update_rank] = column_index;
  }

}// namespace randomgen