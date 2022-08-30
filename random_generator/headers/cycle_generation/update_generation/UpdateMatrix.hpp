#pragma once

#include <iostream>
#include <memory>
#include <span>
#include <vector>

namespace randomgen {

  /**
   * @brief Stores a set of updates, and provides methods to retrieve them. Include metadata about the updates,
   * such as the number of splitting updates, etc..
   *
   *
   * The updates are stored as a matrix of size n_update x update_length.
   */
  class UpdateMatrix {
  private:
    /**
     * @brief Base class for the (Const) iterators. Provides facilities to traverse the update matrix.
     * @tparam UMatrix the type of the update matrix, for const correctness
     */
    template<typename UMatrix>
    class BaseIterator {
    public:
      BaseIterator(UMatrix &parent_matrix, size_t update_index)
          : update_index(update_index), parent_matrix(&parent_matrix) {}

      BaseIterator operator+(int count) const {
        BaseIterator res(*this);
        res += count;
        return res;
      }

      BaseIterator &operator+=(int count) {
        update_index += count;
        return *this;
      }

      BaseIterator &operator++() { return *this += 1; }

      // We can implement the subtraction operators using the addition operators
      BaseIterator operator-(int count) const { return *this + (-count); }


      BaseIterator &operator-=(int count) { return *this += -count; }

      BaseIterator operator--() { return *this -= 1; }

    protected:
      uint32_t update_index;
      UMatrix *parent_matrix;
    };

  public:
    /**
     * @brief Const iterator on the updates
     */
    class ConstIterator : public BaseIterator<const UpdateMatrix> {
    public:
      using BaseIterator::BaseIterator;

      std::span<const double> operator*() const {
        if (update_index >= parent_matrix->n_update or update_index < 0) {
          throw std::out_of_range(
                  "UpdateMatrix::ConstIterator::operator*(): update_index out of range");
        }
        return parent_matrix->getUpdateSpan(update_index);
      }
    };

    /**
     * @brief Const iterator on the updates
     */
    class Iterator : public BaseIterator<UpdateMatrix> {
    public:
      using BaseIterator::BaseIterator;

      std::span<double> operator*() const {
        if (update_index >= parent_matrix->n_update or update_index < 0) {
          throw std::out_of_range(
                  "UpdateMatrix::ConstIterator::operator*(): update_index out of range");
        }
        return parent_matrix->getUpdateSpan(update_index);
      }
    };

    Iterator begin() { return {*this, 0}; }

    ConstIterator begin() const { return {*this, 0}; }

    Iterator end() { return {*this, n_update}; }

    ConstIterator end() const { return {*this, n_update}; }

    UpdateMatrix() = default;

    /**
     * @brief Construct a new matrix able to store updates of the given size
     * @param n_update The total number of updates in the matrix
     * @param n_splits The number of splitting updates
     * @param update_length The length of each update
     */
    UpdateMatrix(uint32_t n_update, uint32_t n_splits, uint32_t update_length);

    UpdateMatrix(const UpdateMatrix &other);
    UpdateMatrix &operator=(const UpdateMatrix &other);

    UpdateMatrix(UpdateMatrix &&other) = default;
    UpdateMatrix &operator=(UpdateMatrix &&other) = default;


    double *getRawUpdates() { return matrix.get(); }
    const double *getRawUpdates() const { return matrix.get(); }

    void setUpdate(uint32_t update_rank, const std::span<const double> &update,
                   uint32_t column_index);

    std::span<double> getUpdateSpan(uint32_t index) {
      std::span<double> res(getRawUpdates() + (index * update_length), update_length);
      return res;
    }

    std::span<const double> getUpdateSpan(uint32_t index) const {
      std::span<const double> res(getRawUpdates() + (index * update_length), update_length);
      return res;
    }

    uint32_t getUpdateCount() const { return n_update; }

    uint32_t getUpdateLength() const { return update_length; }

    uint32_t getSplitCount() const { return n_splits; }

    uint32_t *getUpdateIndices() const { return update_index.get(); }

  private:
    // The number of update in this cycle, and the number of updates that are splits
    // Note that n_update > 0, and n_splits < n_update
    uint32_t n_update = 0, n_splits = 0;
    uint32_t update_length = 0;

    // We use unique ptr to limit the size of this class
    // Thanks to this, this class fits into 32 bytes
    std::unique_ptr<double[]> matrix;
    std::unique_ptr<uint32_t[]> update_index;
  };
}// namespace randomgen
