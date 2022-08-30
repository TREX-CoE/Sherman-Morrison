

#include "UpdateCountExecutor.hpp"

#include <iostream>
#include <omp.h>
#include <vector>

#include "ChainedUpdatePlanner.hpp"
#include "CycleGenerator.hpp"
#include "Engine.hpp"
#include "hdf5CycleOutputStream.hpp"

using namespace randomgen;

size_t dumpCyclesToDatastream(hdf5CycleOutputStream &ds, std::vector<Cycle> &cycle_buffer) {
#pragma omp critical
  for (auto &c: cycle_buffer) { ds << c; }

  // Number of cycles that were dumped
  size_t res = cycle_buffer.size();
  cycle_buffer.clear();
  return res;
}


void generateUpdates(std::filesystem::path output_path, const size_t matrix_size,
                     const size_t y_padding, const size_t lower_update_count,
                     const size_t upper_update_count, const size_t max_split_count,
                     const size_t cycle_per_rank) {


  auto planner = std::make_shared<ChainedUpdatePlanBuilder>();
  auto matrix_generator = std::make_shared<UpdateMatrixGenerator>(1e-6, 1e-6);
  auto update_generator = std::make_shared<UpdateGenerator>(planner, matrix_generator);
  CycleGenerator generator(matrix_size, y_padding, update_generator);

  hdf5CycleOutputStream ds(output_path);

  // For pretty printing the progress
  // The rank is a combination of (n_update, n_splits)
  size_t total_rank = 0;
  size_t current_rank = 0;
  size_t cycle_counter = 0;

  // Total number of ranks to compute
  for (size_t i = lower_update_count; i <= upper_update_count; i++) {
    size_t splits_count = std::min(max_split_count, i);
    splits_count = std::max(splits_count - 2, 0UL);
    total_rank += splits_count;
  }


#pragma omp parallel default(none)                                                                 \
        shared(total_rank, current_rank, cycle_counter, generator, ds, std::cout,                  \
               lower_update_count, upper_update_count, max_split_count, cycle_per_rank)
  {
    // Generate this amount of cycles before dumping them
    // This is to avoid contention
    const size_t cycle_buffer_size = std::min(256UL, 2 * cycle_per_rank);
    std::vector<Cycle> cycle_buffer;

    // Some cycles take longer to generate than others
    // (espcially cycles with more splits and / or more updates)
    // So we use a dynamic scheduling
#pragma omp for schedule(dynamic)
    for (size_t n_update = lower_update_count; n_update <= upper_update_count; n_update++) {
      size_t splits_count = std::min(max_split_count, n_update - 1);
      // For each update count, generate all combinations of splits number
      for (size_t n_splits = 2; n_splits <= splits_count; n_splits++) {

        // Pretty print the progress :)
        if (omp_get_thread_num() == 0) {
          std::cout << "Generated rank " << current_rank << "/" << total_rank << " ("
                    << cycle_counter << ")" << std::endl;
        }

        for (size_t k = 0; k < cycle_per_rank; k++) {
          cycle_buffer.push_back(generator.make(n_update, n_splits));
          // If we have generated enough cycles, dump them

          if (cycle_buffer.size() == cycle_buffer_size) {
            // ensure the counter increase is atomic
            size_t tmp = dumpCyclesToDatastream(ds, cycle_buffer);
#pragma omp atomic
            cycle_counter += tmp;
          }
        }
#pragma omp atomic
        current_rank++;
      }
    }
    // Dump the remaining cycles still in the buffer
    size_t tmp = dumpCyclesToDatastream(ds, cycle_buffer);
#pragma omp atomic
    cycle_counter += tmp;
  }
  std::cout << "Generated rank " << current_rank << "/" << total_rank << " (" << cycle_counter
            << ")" << std::endl;
}


int UpdateCountExecutor::exec() {

  // The generated matrix are square matrices...
  // But we can pad them with zeros to make them rectangular (useful for optimal vectorization)
  constexpr size_t matrix_size = 21;
  constexpr size_t y_padding = 3;

  // Start generating cycles with this number of updates
  constexpr size_t lower_update_count = 1;
  // The maximal number of updates per cycle to generate
  constexpr size_t upper_update_count = 20;
  // The maximum number of splits per cycle
  // By default, this is capped to upper_update_count - 1
  constexpr size_t max_split_count = upper_update_count - 1;
  // The number of cycles to generate for each "rank"
  // Where the rank is the combination (update_count, split_count)
  constexpr size_t cycle_per_rank = 80;

  auto output_path = parent->getOutputPath();

  generateUpdates(output_path, matrix_size, y_padding, lower_update_count, upper_update_count,
                  max_split_count, cycle_per_rank);
  return 0;
}
