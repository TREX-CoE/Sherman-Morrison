#include <iostream>
#include <string>
#include "hdf5/serial/hdf5.h"
#include "hdf5/serial/H5Cpp.h"

#include "SM_MaponiA3.hpp"
#include "SM_Standard.hpp"
#include "Helpers.hpp"

using namespace H5;
//#define DEBUG

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

  unsigned int dim, nupdates, col, i, j;
  read_int(file, group + "/slater_matrix_dim", &dim);
  read_int(file, group + "/nupdates", &nupdates);

  double * slater_matrix = new double[dim*dim];
  read_double(file, group + "/slater_matrix", slater_matrix);

  double * slater_inverse = new double[dim*dim];
  read_double(file, group + "/slater_inverse", slater_inverse);
  slater_inverse = transpose(slater_inverse, dim);

  unsigned int * col_update_index = new unsigned int[nupdates];
  read_int(file, group + "/col_update_index", col_update_index);

  double * updates = new double[nupdates*dim];
  read_double(file, group + "/updates", updates);

  /* Test */
#ifdef DEBUG
  showMatrix(slater_matrix, dim, "OLD Slater");
#endif

#ifdef DEBUG
  showMatrix(slater_inverse, dim, "OLD Inverse");
#endif

  for (j = 0; j < nupdates; j++) {
    for (i = 0; i < dim; i++) {
      col = col_update_index[j];
      slater_matrix[i*dim + (col - 1)] += updates[i + j*dim];
    }
  }

  MaponiA3(slater_inverse, dim, nupdates, updates, col_update_index);
  //SM(slater_inverse, dim, nupdates, updates, col_update_index);

#ifdef DEBUG
  showMatrix(slater_matrix, dim, "NEW Slater");
#endif

#ifdef DEBUG
  showMatrix(slater_inverse, dim, "NEW Inverse");
#endif

  double * res = new double[dim*dim] {0};
  matMul(slater_matrix, slater_inverse, res, dim);
  bool ok = is_identity(res, dim, 0.5e-4);

#ifdef DEBUG
  showMatrix(res, dim, "Result");
#endif

  delete [] res, updates, col_update_index,
            slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Execute from within 'datasets/'" << std::endl;
    std::cerr << "usage: test_internal_h5 <start cycle> <stop cycle>" << std::endl;
    return 1;
  }
  int start_cycle = std::stoi(argv[1]);
  int stop_cycle = std::stoi(argv[2]);
  H5File file(FILE_NAME, H5F_ACC_RDONLY);

  bool ok;
  for (int cycle = start_cycle; cycle < stop_cycle+1; cycle++) {
    ok = test_cycle(file, cycle);
    if (ok) {
      std::cerr << "ok -- cycle " << std::to_string(cycle)
      << std::endl;
    }
    else {
      std::cerr << "failed -- cycle " << std::to_string(cycle)
      << std::endl;
    }
  }

  return ok;
}
