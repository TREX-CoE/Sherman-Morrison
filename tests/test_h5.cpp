#include "hdf5/serial/H5Cpp.h"
#include "hdf5/serial/hdf5.h"

#include "Helpers.hpp"
#include "SM_Maponi.hpp"
#include "SM_Standard.hpp"
#include "SMWB.hpp"

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
#ifdef DEBUG
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

#ifdef DEBUG
  showMatrix(slater_matrix, dim, "OLD Slater");
#endif

#ifdef DEBUG
  showMatrix(u, dim, "Updates");
#endif

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
  } else if (version == "smwb1") {
    SMWB1(slater_inverse, dim, nupdates, u, col_update_index);
#ifdef MKL
  } else if (version == "lapack") {
    memcpy(slater_inverse, slater_matrix, dim * dim * sizeof(double));
    inverse(slater_inverse, dim);
#endif
  } else {
    std::cerr << "Unknown version " << version << std::endl;
    exit(1);
  }

#ifdef DEBUG
  showMatrix(slater_matrix, dim, "NEW Slater");
#endif

#ifdef DEBUG
  showMatrix(slater_inverse, dim, "NEW Inverse");
#endif

  double *res = new double[dim * dim]{0};
  matMul(slater_matrix, slater_inverse, res, dim);
  bool ok = is_identity(res, dim, tolerance);

  double res_max = residual_max(res, dim);
  double res2 = residual_frobenius2(res, dim);
  // double det;
  // double **tmp = new double *[dim];
  // for (int i = 0; i < dim; i++) {
  //   tmp[i] = new double[dim];
  //   for (int j = 0; j < dim; j++) {
  //     tmp[i][j] = res[i * dim + j];
  //   }
  // }
  // det = determinant(tmp, dim);
  // delete[] tmp;
  // std::cout << "Residual = " << version << " " << cycle << " " << res_max <<
  // " "
  //           << res2 << " " << det << std::endl;
  std::cout << "Residual = " << version << " " << cycle << " " << res_max << " "
            << res2 << std::endl;

#ifdef DEBUG
  showMatrix(res, dim, "Result");
#endif

  delete[] res, updates, u, col_update_index, slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr
        << "usage: test_h5 <version> <start cycle> <stop cycle> <tolerance>"
        << std::endl;
    return 1;
  }
  std::string version(argv[1]);
  int start_cycle = std::stoi(argv[2]);
  int stop_cycle = std::stoi(argv[3]);
  double tolerance = std::stod(argv[4]);
  H5File file(FILE_NAME, H5F_ACC_RDONLY);

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
