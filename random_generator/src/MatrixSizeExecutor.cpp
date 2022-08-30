
#include "MatrixSizeExecutor.hpp"

#include <iostream>
#include <omp.h>
#include <vector>

#include "ChainedUpdatePlanner.hpp"
#include "CycleGenerator.hpp"
#include "SizedBasedCycleStream.hpp"

#include "Engine.hpp"

using namespace randomgen;

namespace {

  size_t dumpCyclesToDatastream(hdf5CycleOutputStream &ds, std::vector<Cycle> &cycle_buffer) {
#pragma omp critical
    for (auto &c: cycle_buffer) { ds << c; }

    // Number of cycles that were dumped
    size_t res = cycle_buffer.size();
    cycle_buffer.clear();
    return res;
  }


  void checkMultisizedParameters(const size_t starting_matrix_size, const size_t final_matrix_size,
                                 const size_t update_count, const size_t max_split_count,
                                 const size_t cycle_per_rank) {
    if (starting_matrix_size > final_matrix_size) {
      throw std::runtime_error("Starting matrix size must be smaller than final matrix size");
    }

    if (max_split_count > update_count) {
      throw std::runtime_error("Update count must be smaller than max split count");
    }

    if (cycle_per_rank == 0) { throw std::runtime_error("Cycle per rank must be greater than 0"); }

    if (update_count > starting_matrix_size) {
      throw std::runtime_error("Update count must be smaller than starting matrix size");
    }
  }

  void generate_multisized_dataset(const std::filesystem::path &output_path,
                                   const size_t starting_matrix_size,
                                   const size_t final_matrix_size, const size_t update_count,
                                   const size_t max_split_count, const size_t cycle_per_rank) {
    checkMultisizedParameters(starting_matrix_size, final_matrix_size, update_count,
                              max_split_count, cycle_per_rank);


    auto planner = std::make_shared<ChainedUpdatePlanBuilder>();
    auto matrix_generator = std::make_shared<UpdateMatrixGenerator>(1e-3, 1e-6);
    auto update_generator = std::make_shared<UpdateGenerator>(planner, matrix_generator);

    SizedBasedCycleStream ds(output_path);

    size_t total_rank_count = 0;
    for (size_t i = starting_matrix_size; i <= final_matrix_size; i++) {
      long splits_count = std::min(max_split_count, update_count - 1);
      splits_count = std::max(splits_count, 0L);
      total_rank_count += (splits_count + 1);
    }
    size_t generated_rank = 0;
    size_t generated_cycles = 0;


#pragma omp parallel default(none) shared(generated_cycles, generated_rank)                        \
        shared(starting_matrix_size, final_matrix_size, update_generator, cycle_per_rank,          \
               max_split_count, update_count, ds, total_rank_count, std::cout)
    {
      // Generate this amount of cycles before dumping them
      // This is to avoid contention
      const size_t cycle_buffer_size = std::min(256UL, 2 * cycle_per_rank);
      std::vector<Cycle> cycle_buffer;


      // Iterate on the matrix size
#pragma omp for schedule(dynamic)
      for (size_t curr_size = starting_matrix_size; curr_size <= final_matrix_size; curr_size++) {
        CycleGenerator generator(curr_size, 0, update_generator);
        size_t splits_count = std::min(max_split_count, update_count - 1);

        // Iterate on the number of splits
        for (size_t n_splits = 0; n_splits <= splits_count; n_splits++) {

          if (omp_get_thread_num() == 0) {
            std::cout << "Generated " << generated_rank << "/" << total_rank_count << " ranks"
                      << "(" << generated_cycles << ")" << std::endl;
          }

          // Generate cycle_per_rank matrices
          for (size_t k = 0; k < cycle_per_rank; k++) {
            cycle_buffer.push_back(generator.make(update_count, n_splits));
            // If we have generated enough cycles, dump them

            if (cycle_buffer.size() == cycle_buffer_size) {
              // ensure the counter increase is atomic
              size_t tmp = dumpCyclesToDatastream(ds, cycle_buffer);
#pragma omp atomic
              generated_cycles += tmp;
            }
          }
#pragma omp atomic
          generated_rank += 1;
        }
      }
    }
  }
}// namespace

int MatrixSizeExecutor::exec() {

  constexpr size_t starting_matrix_size = 31;
  constexpr size_t ending_matrix_size = 200;
  constexpr size_t nbr_update = 31;
  constexpr size_t max_split_count = nbr_update - 1;
  constexpr size_t cycle_per_rank = 20;

  auto output_path = parent->getOutputPath();

  generate_multisized_dataset(output_path, starting_matrix_size, ending_matrix_size, nbr_update,
                              max_split_count, cycle_per_rank);

  return 0;
}