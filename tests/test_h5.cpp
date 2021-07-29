#include "hdf5/serial/H5Cpp.h"
#include "hdf5/serial/hdf5.h"

#include "Helpers.hpp"
#include "SMWB.hpp"
#include "SM_Maponi.hpp"
#include "SM_Standard.hpp"
#include "Woodbury.hpp"

#define PERF

#ifdef PERF
unsigned int repetition_number;
#endif

const H5std_string FILE_NAME("dataset.hdf5");

void read_int(H5::H5File file, std::string key, unsigned int *data) {
  H5::DataSet ds = file.openDataSet(key);
  ds.read(data, H5::PredType::STD_U32LE);
  ds.close();
}

void read_double(H5::H5File file, std::string key, double *data) {
  H5::DataSet ds = file.openDataSet(key);
  ds.read(data, H5::PredType::IEEE_F64LE);
  ds.close();
}

int test_cycle(H5::H5File file, int cycle, std::string version, double breakdown, double tolerance) {

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

  // Transform replacement updates in 'updates[]' into additive updates in 'u[]'
  for (j = 0; j < nupdates; j++) {
    for (i = 0; i < dim; i++) {
      col = col_update_index[j];
      u[i + j * dim] =
          updates[i + j * dim] - slater_matrix[i * dim + (col - 1)];
      slater_matrix[i * dim + (col - 1)] = updates[i + j * dim];
    }
  }
  delete[] updates;

#ifdef PERF
  std::cout << "# of reps. = " << repetition_number << std::endl;
  double *slater_inverse_nonpersistent = new double[dim * dim];

  if (version == "sm1") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      SM1(slater_inverse_nonpersistent, dim, nupdates,
        u, col_update_index, breakdown);
    }
  }
  else if (version == "wb2") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      WB2(slater_inverse_nonpersistent, dim,
        u, col_update_index, breakdown);
    }
  }
  else if (version == "wb3") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      WB3(slater_inverse_nonpersistent, dim,
        u, col_update_index, breakdown);
    }
  }
  else if (version == "sm2") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      SM2(slater_inverse_nonpersistent, dim, nupdates,
        u, col_update_index, breakdown);
    }
  }
  else if (version == "wb2s") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      WB2s(slater_inverse_nonpersistent, dim, nupdates,
        u, col_update_index, breakdown);
    }
  }
  else if (version == "wb3s") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      WB3s(slater_inverse_nonpersistent, dim, nupdates,
        u, col_update_index, breakdown);
    }
  }
  else if (version == "wb32s") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      WB32s(slater_inverse_nonpersistent, dim, nupdates,
        u, col_update_index, breakdown);
    }
  }
#ifdef MKL
  else if (version == "lapack") {
    memcpy(slater_inverse_nonpersistent, slater_matrix,
            dim * dim * sizeof(double));
    inverse(slater_inverse_nonpersistent, dim);
  }
#endif // MKL
  else {
    std::cerr << "Unknown version " << version << std::endl;
    exit(1);
  }
  std::memcpy(slater_inverse, slater_inverse_nonpersistent,
              dim * dim * sizeof(double));
  delete[] slater_inverse_nonpersistent;
#else
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
  } else if (version == "wb2s") {
    WB2s(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "wb3s") {
    WB3s(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "wb32s") {
    WB32s(slater_inverse, dim, nupdates, u, col_update_index);
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
  delete[] u, col_update_index;

  showMatrix(slater_matrix, dim, "Slater Matrix");
  showMatrix(slater_inverse, dim, "Slater Inverse");
  double *res = new double[dim * dim]{0};
  {
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      for (unsigned int k = 0; k < dim; k++) {
        res[i * dim + j] += slater_matrix[i * dim + k] * slater_inverse[k * dim + j];
      }
    }
  }
  }

  //matMul2(slater_matrix, slater_inverse, res, dim, dim, dim);
  //
  //
  for (unsigned int i = 0; i < dim; i++) {
    printf("[");
    for (unsigned int j = 0; j < dim; j++) {
      if (slater_matrix[i * dim + j] >= 0) {
        printf("  %17.10e,", slater_matrix[i * dim + j]);
      } else {
        printf(" %17.10e,", slater_matrix[i * dim + j]);
      }
    }
    printf(" ],\n");
  }
  printf("\n\n");
  //
  //
  //
  //
  for (unsigned int i = 0; i < dim; i++) {
    printf("[");
    for (unsigned int j = 0; j < dim; j++) {
      if (slater_inverse[i * dim + j] >= 0) {
        printf("  %17.10e,", slater_inverse[i * dim + j]);
      } else {
        printf(" %17.10e,", slater_inverse[i * dim + j]);
      }
    }
    printf(" ],\n");
  }
  printf("\n\n");
  //
  //
  //
  //
  for (unsigned int i = 0; i < dim; i++) {
    printf("[");
    for (unsigned int j = 0; j < dim; j++) {
      if (res[i * dim + j] >= 0) {
        printf("  %17.10e,", res[i * dim + j]);
      } else {
        printf(" %17.10e,", res[i * dim + j]);
      }
    }
    printf(" ],\n");
  }
  printf("\n\n");
  //
  //
  bool ok = is_identity(res, dim, tolerance);
  double res_max = residual_max(res, dim);
  double res2 = residual_frobenius2(res, dim);

  std::cout << "Residual = " << version << " " << cycle << " " << res_max << " "
            << res2 << std::endl;

#ifdef DEBUG2
  showMatrix(res, dim, "Result");
#endif

  delete[] res, slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
#ifdef PERF
  if (argc != 7) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <start cycle> <stop cycle> <break-down threshold> <tolerance> <number of reps.>"
        << std::endl;
    return 1;
  }
#else
  if (argc != 6) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <start cycle> <stop cycle> <break-down threshold> <tolerance>"
        << std::endl;
    return 1;
  }
#endif

  std::string version(argv[1]);
  int start_cycle = std::stoi(argv[2]);
  int stop_cycle = std::stoi(argv[3]);
  double breakdown = std::stod(argv[4]);
  double tolerance = std::stod(argv[5]);
  H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

#ifdef PERF
  repetition_number = std::stoi(argv[6]);
#endif

  bool ok;
  for (int cycle = start_cycle; cycle < stop_cycle + 1; cycle++) {
    ok = test_cycle(file, cycle, version, breakdown, tolerance);
    if (ok) {
      std::cerr << "ok -- cycle " << std::to_string(cycle) << std::endl;
    } else {
      std::cerr << "failed -- cycle " << std::to_string(cycle) << std::endl;
    }
  }

  return ok;
}
