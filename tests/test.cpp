#include <iostream>
#include <string>
#include "hdf5/serial/hdf5.h"
#include "H5Cpp.h"

#include "../SM_MaponiA3.hpp"
#include "../Helpers.hpp"

using namespace H5;
#define DEBUG 1

const H5std_string FILE_NAME( "datasets.hdf5" );

void read_int(H5File file, std::string key, unsigned int * data) {
  DataSet ds = file.openDataSet(key);
  ds.read(data, PredType::STD_U32LE);
  ds.close();
}

void read_double(H5File file, std::string key, double * data) {
  DataSet ds = file.openDataSet(key);
  ds.read(data, PredType::IEEE_F64LE);
  ds.close();
}

int test_cycle(H5File file, int cycle) {

  /* Read the data */

  std::string group = "cycle_" + std::to_string(cycle);

  unsigned int dim, nupdates;
  read_int(file, group + "/slater_matrix_dim", &dim);
  read_int(file, group + "/nupdates", &nupdates);

  double * slater_matrix = new double[dim*dim];
  read_double(file, group + "/slater_matrix", slater_matrix);

  double * slater_inverse = new double[dim*dim];
  read_double(file, group + "/slater_inverse", slater_inverse);

  unsigned int * col_update_index = new unsigned int[nupdates];
  read_int(file, group + "/col_update_index", col_update_index);

  double * updates = new double[nupdates*dim];
  read_double(file, group + "/updates", updates);

  /* Test */
#ifdef DEBUG
  showMatrix(slater_matrix, dim, "Slater");
#endif

  MaponiA3(slater_inverse, dim, nupdates, updates, col_update_index);

#ifdef DEBUG
  showMatrix(slater_inverse, dim, "Inverse");
#endif

  double * res = matMul(slater_matrix, slater_inverse, dim);
  bool ok = is_identity(res, dim, 1.0e-8);

#ifdef DEBUG
  showMatrix(res, dim, "Result");
#endif

  delete [] res, updates, col_update_index, slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: ./test <cycle>" << std::endl;
    return 1;
  }
  int cycle = std::stoi(argv[1]);
  H5File file(FILE_NAME, H5F_ACC_RDONLY);

  bool ok = test_cycle(file, cycle);

  if (ok) {
    std::cerr << "ok -- cycle " << std::to_string(cycle) << std::endl;
  } else {
    std::cerr << "failed -- cycle " << std::to_string(cycle) << std::endl;
  }
  return ok;
}

