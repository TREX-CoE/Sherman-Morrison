#include "hdf5/serial/H5Cpp.h"
#include "hdf5/serial/hdf5.h"

#include "Helpers.hpp"
#include "SM_Maponi.hpp"
#include "SM_Standard.hpp"
#include "Woodbury.hpp"
#include "SMWB.hpp"
#include <fstream>
#include <vector>

#define PERF
// #define STATUS
// #define RESIDUAL

#ifdef PERF
unsigned int repetition_number;
#endif

using namespace H5;

// #define DEBUG

const H5std_string FILE_NAME("dataset.hdf5");

void read_int(H5File file, std::string key, unsigned int *data) {
  DataSet ds = file.openDataSet(key);
  ds.read(data, PredType::STD_U32LE);
  ds.close();
}

void read_double(H5File file, std::string key, double *data) {
  DataSet ds = file.openDataSet(key);
  ds.read(data, PredType::IEEE_F64LE);
  ds.close();
}

int test_cycle(H5File file, int cycle, std::string version, double tolerance) {

  /* Read the data */

  std::string group = "cycle_" + std::to_string(cycle);

  unsigned int dim, nupdates, col, i, j;
  read_int(file, group + "/slater_matrix_dim", &dim);
  read_int(file, group + "/nupdates", &nupdates);

  double *slater_matrix = new double[dim * dim];
  read_double(file, group + "/slater_matrix", slater_matrix);

  double *slater_inverse = new double[dim * dim];
  read_double(file, group + "/slater_inverse", slater_inverse);

  unsigned int *col_update_index = new unsigned int[nupdates];
  read_int(file, group + "/col_update_index", col_update_index);

  double *updates = new double[nupdates * dim];
  read_double(file, group + "/updates", updates);

  double *u = new double[nupdates * dim];

  /* Test */
  #ifdef DEBUG2
  showMatrix(slater_inverse, dim, "OLD Inverse");
  #endif

  // Transform replacement updates in 'updates[]' into additive updates in 'u[]'
  for (j = 0; j < nupdates; j++) {
    for (i = 0; i < dim; i++) {
      col = col_update_index[j];
      u[i + j * dim] =
          updates[i + j * dim] - slater_matrix[i * dim + (col - 1)];
      slater_matrix[i * dim + (col - 1)] = updates[i + j * dim];
    }
  }

  #ifdef DEBUG2
  showMatrix(slater_matrix, dim, "OLD Slater");
  showMatrix(u, dim, "Updates");
  #endif

  #ifdef PERF
  #ifdef DEBUG1
  std::cerr << "# of reps. = " << repetition_number << std::endl;
  #endif // DEBUG1
  double *slater_inverse_nonpersistent = new double[dim * dim];
  if (version == "sm1") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      SM1(slater_inverse_nonpersistent, dim, nupdates, u, col_update_index);
    }
  } else if (version == "sm2") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      SM2(slater_inverse_nonpersistent, dim, nupdates, u, col_update_index);
    }
  } else if (version == "sm3") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      SM3(slater_inverse_nonpersistent, dim, nupdates, u, col_update_index);
    }
  } else if (version == "sm4") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      SM4(slater_inverse_nonpersistent, dim, nupdates, u, col_update_index);
    }
  } else if (version == "wb2") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      WB2(slater_inverse_nonpersistent, dim, u, col_update_index);
    }
  } else if (version == "wb3") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      WB3(slater_inverse_nonpersistent, dim, u, col_update_index);
    }
  } else if (version == "smwb1") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      SMWB1(slater_inverse_nonpersistent, dim, nupdates, u, col_update_index);
    }
  } else if (version == "smwb4") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_inverse, dim * dim * sizeof(double));
      SMWB4(slater_inverse_nonpersistent, dim, nupdates, u, col_update_index);
    }
  #ifdef MKL
  } else if (version == "lapack") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      std::memcpy(slater_inverse_nonpersistent, slater_matrix, dim * dim * sizeof(double));
      inverse(slater_inverse_nonpersistent, dim);
    }
  #endif // MKL
  } else {
    std::cerr << "Unknown version " << version << std::endl;
    exit(1);
  }
  std::memcpy(slater_inverse, slater_inverse_nonpersistent, dim * dim * sizeof(double));
  delete[] slater_inverse_nonpersistent;
  #else //  No performance measurements repetition
  if (version == "maponia3") {
    MaponiA3(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "maponia3s") {
    MaponiA3S(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm1") {
    SM1(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm2") {
    SM2(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm3") {
    SM3(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm4") {
    SM4(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "wb2") {
    WB2(slater_inverse, dim, u, col_update_index);
  } else if (version == "wb3") {
    WB3(slater_inverse, dim, u, col_update_index);
  } else if (version == "smwb1") {
    SMWB1(slater_inverse, dim, nupdates, u, col_update_index);
  // } else if (version == "smwb2") {
  //   SMWB2(slater_inverse, dim, nupdates, u, col_update_index);
  // } else if (version == "smwb3") {
  //   SMWB3(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "smwb4") {
    SMWB4(slater_inverse, dim, nupdates, u, col_update_index);
  #ifdef MKL
  } else if (version == "lapack") {
    memcpy(slater_inverse, slater_matrix, dim * dim * sizeof(double));
    inverse(slater_inverse, dim);
  #endif // MKL
  } else {
    std::cerr << "Unknown version " << version << std::endl;
    exit(1);
  }
  #endif // PERF

  #ifdef DEBUG2
  showMatrix(slater_matrix, dim, "NEW Slater");
  showMatrix(slater_inverse, dim, "NEW Inverse");
  #endif

  double *res = new double[dim * dim]{0};
  matMul(slater_matrix, slater_inverse, res, dim);
  bool ok = is_identity(res, dim, tolerance);
  double res_max = residual_max(res, dim);
  double res2 = residual_frobenius2(res, dim);

  #ifdef RESIDUAL
  std::cout << "Residual = " << version << " " << cycle << " " << res_max << " "
            << res2 << std::endl;
  #endif

  #ifdef DEBUG2
  showMatrix(res, dim, "Result");
  #endif

  delete[] res, updates, u, col_update_index, slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
  #ifdef PERF
  if (argc != 5) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <cycle file> <tolerance> <number of reps.>"
        << std::endl;
    return 1;
  }
  #else
  if (argc != 4) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <cycle file> <tolerance>"
        << std::endl;
    return 1;
  }
  #endif    
  std::string version(argv[1]);
  std::string cyclefile_name(argv[2]);
  std::ifstream cyclefile(cyclefile_name);
  std::vector<int> cycles;
  unsigned int cycle;
  while (cyclefile >> cycle) cycles.push_back(cycle);
  double tolerance = std::stod(argv[3]);
  H5File file(FILE_NAME, H5F_ACC_RDONLY);

  #ifdef PERF
  repetition_number = std::stoi(argv[4]);
  #endif

  bool ok;
  for (auto & cycle : cycles) {
    ok = test_cycle(file, cycle, version, tolerance);
    #ifdef STATUS
    if (ok) {
      std::cerr << "ok -- cycle " << std::to_string(cycle) << std::endl;
    } else {
      std::cerr << "failed -- cycle " << std::to_string(cycle) << std::endl;
    }
    #endif
  }

  return ok;
}
