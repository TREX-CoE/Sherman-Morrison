#include "hdf5/serial/H5Cpp.h"
#include "hdf5/serial/hdf5.h"

#include "Helpers.hpp"

extern "C" {
#include "qmckl.h"
}

#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#ifndef THRESHOLD
#define THRESHOLD 1e-3
#endif

#include "cstring"
#include "iostream"

#define PERF

#ifdef PERF
unsigned int repetition_number;
#endif

const H5std_string FILE_NAME("dataset.hdf5");

void read_int(H5::H5File file, std::string key, uint64_t *data) {
  H5::DataSet ds = file.openDataSet(key);
  ds.read(data, H5::PredType::STD_U32LE);
  ds.close();
}

void read_double(H5::H5File file, std::string key, double *data) {
  H5::DataSet ds = file.openDataSet(key);
  ds.read(data, H5::PredType::IEEE_F64LE);
  ds.close();
}

int test_cycle(H5::H5File file, int cycle, std::string version, double tolerance) {

  /* Read the data */

  std::string group = "cycle_" + std::to_string(cycle);

  unsigned int col, i, j;
  uint64_t dim, nupdates;

  read_int(file, group + "/slater_matrix_dim", &dim);
  read_int(file, group + "/nupdates", &nupdates);

  double *slater_matrix = new double[dim * dim];
  read_double(file, group + "/slater_matrix", slater_matrix);

  double *slater_inverse = new double[dim * dim];
  read_double(file, group + "/slater_inverse", slater_inverse);

  uint64_t *col_update_index = new uint64_t[nupdates];
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
  double *slater_inverse_nonpersistent = new double[dim * dim];
  if (version == "sm1") {
    for (unsigned int i = 0; i < repetition_number; i++) {
      memcpy(slater_inverse_nonpersistent, slater_inverse,
                  dim * dim * sizeof(double));
      qmckl_context context;
      context = qmckl_context_create();
      qmckl_exit_code rc;
      // const uint64_t cdim = dim;
      // const uint64_t cnupdates = nupdates;
      // const uint64_t* ccol_update_index = col_update_index;
      // const double* cu = u;
      // rc = qmckl_sherman_morrison_c(
      //   context,
      //   cdim,
      //   cnupdates,
      //   cu,
      //   ccol_update_index,
      //   slater_inverse_nonpersistent);
      const uint64_t Dim = 2;
      const uint64_t N_updates = 2;
      const uint64_t Updates_index[2] = {0, 0};
      const double Updates[4] = {0.0, 0.0, 0.0, 0.0};
      double Slater_inv[4] = {0.0, 0.0, 0.0, 0.0};
      rc = qmckl_sherman_morrison_c(context,
            Dim,
            N_updates,
            Updates,
            Updates_index,
            Slater_inv);
      assert(rc == QMCKL_SUCCESS);
    }
  }
  else {
    std::cerr << "Unknown version " << version << std::endl;
    exit(1);
  }
  std::memcpy(slater_inverse, slater_inverse_nonpersistent,
              dim * dim * sizeof(double));
  delete[] slater_inverse_nonpersistent;
#else //  No performance measurements repetition
  if (version == "sm1") {
      qmckl_exit_code sherman_morrison_exit;
      qmckl_context context;
      sherman_morrison_exit = qmckl_sherman_morrison_c(context,
        dim,
        nupdates,
        u,
        col_update_index,
        slater_inverse);
  }
  else {
    std::cerr << "Unknown version " << version << std::endl;
    exit(1);
  }
#endif // PERF
  delete[] u, col_update_index;

  double *res = new double[dim * dim]{0};
  matMul(slater_matrix, slater_inverse, res, dim);
  bool ok = is_identity(res, dim, tolerance);
  double res_max = residual_max(res, dim);
  double res2 = residual_frobenius2(res, dim);

  std::cout << "Residual = " << version << " " << cycle << " " << res_max << " "
            << res2 << std::endl;

  delete[] res, slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
#ifdef PERF
  if (argc != 6) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <start cycle> <stop cycle> <tolerance> <number of reps.>"
        << std::endl;
    return 1;
  }
#else
  if (argc != 5) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <start cycle> <stop cycle> <tolerance>"
        << std::endl;
    return 1;
  }
#endif

  std::string version(argv[1]);
  int start_cycle = std::stoi(argv[2]);
  int stop_cycle = std::stoi(argv[3]);
  double tolerance = std::stod(argv[4]);
  H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

#ifdef PERF
  repetition_number = std::stoi(argv[5]);
#endif

  bool ok;
  for (int cycle = start_cycle; cycle < stop_cycle + 1; cycle++) {
    ok = test_cycle(file, cycle, version, tolerance);
    if (ok) {
      std::cerr << "ok -- cycle " << std::to_string(cycle) << std::endl;
    } else {
      std::cerr << "failed -- cycle " << std::to_string(cycle) << std::endl;
    }
  }

  return ok;
}
