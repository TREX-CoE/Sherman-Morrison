// This files is almost the same as test_h5.cpp, with the difference that it
// dumps Verificarlo probes for vfc_ci integration, and that it reads a list of
// cycles in a CSV file, instead of accepting a start and an end cycle (which
// makes it easier to select the exact cycles we are interested in with vfc_ci).

#include <hdf5/serial/hdf5.h>
#include <hdf5/serial/H5Cpp.h>

#include <vector>
#include <fstream>
#include <sstream>


#include "SM_MaponiA3.hpp"
#include "SM_Standard.hpp"
#include "SM_Helpers.hpp"
#include "vfc_probe.h"

using namespace H5;
// #define DEBUG

const H5std_string FILE_NAME( "datasets/ci_dataset.hdf5" );

double residual_max(double * A, unsigned int Dim) {
  double max = 0.0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      double delta = (A[i * Dim + j] - (i == j));
      delta = abs(delta);
      if (delta > max) max = delta;
    }
  }
  return max;
}

double residual2(double * A, unsigned int Dim) {
  double res = 0.0;
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      double delta = (A[i * Dim + j] - (i == j));
      res += delta*delta;
    }
  }
  return res;
}

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


/* Return a vector containing all cycles to execute by reading a data file */
std::vector<int> get_cycles_list(std::string path) {
  std::ifstream file_stream(path);
  std::stringstream string_stream;
  string_stream << file_stream.rdbuf();

  std::string cycle_str;
  std::vector<int> cycles_list = {};

  while(string_stream >> cycle_str) {
    cycles_list.push_back(std::stoi(cycle_str));
  }

  return cycles_list;
}

int test_cycle(H5File file, int cycle, std::string version, vfc_probes * probes) {

  /* Read the data */

  std::string group = "cycle_" + std::to_string(cycle);

  // This will result in the same string as group but with the cycle number
  // being zero-padded. This is used when calling vfc_put_probe later on.
  std::string zero_padded_group = std::to_string(cycle);
  zero_padded_group = "cycle_" +
  std::string(5 - zero_padded_group.length(), '0') + zero_padded_group;

  try{
    file.openGroup(group);
  } catch(H5::Exception& e){
    std::cerr << "group " << group << "not found" << std::endl;
    return 0;
  }

  unsigned int dim, nupdates, col, i, j;
  read_int(file, group + "/slater_matrix_dim", &dim);
  read_int(file, group + "/nupdates", &nupdates);


  double * slater_matrix = new double[dim*dim];
  read_double(file, group + "/slater_matrix", slater_matrix);

  double * slater_inverse = new double[dim*dim];
  read_double(file, group + "/slater_inverse", slater_inverse);
  //slater_inverse = transpose(slater_inverse, dim);

  unsigned int * col_update_index = new unsigned int[nupdates];
  read_int(file, group + "/col_update_index", col_update_index);

  double * updates = new double[nupdates*dim];
  read_double(file, group + "/updates", updates);

  double * u = new double[nupdates*dim];

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
      u[i + j*dim] = updates[i + j*dim] - slater_matrix[i*dim + (col - 1)];
      slater_matrix[i*dim + (col - 1)] = updates[i + j*dim];
    }
  }

  if (version == "maponia3") {
    MaponiA3(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm1") {
    SM1(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm2") {
    SM2(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm3") {
    SM3(slater_inverse, dim, nupdates, u, col_update_index);
  } else if (version == "sm3") {
    SM4(slater_inverse, dim, nupdates, u, col_update_index);
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

  double * res = new double[dim*dim] {0};
  matMul(slater_matrix, slater_inverse, res, dim);
  bool ok = is_identity(res, dim, 1e-3);

  double res_max = residual_max(res, dim);
  double res2 = residual2(res, dim);

#ifdef DEBUG
  showMatrix(res, dim, "Result");
#endif

  vfc_put_probe(probes, &(zero_padded_group)[0], &("res_max_" + version)[0], res_max);
  vfc_put_probe(probes, &(zero_padded_group)[0], &("res2_" + version)[0], res2);

  delete [] res, updates, u, col_update_index,
            slater_matrix, slater_inverse;

  return ok;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Execute from within '/'" << std::endl;
    std::cerr << "usage: test_h5 <version> <path to cycles file>" << std::endl;
    return 1;
  }
  std::string version(argv[1]);
  std::vector<int> cycles_list = get_cycles_list(argv[2]);
  H5File file(FILE_NAME, H5F_ACC_RDONLY);

  vfc_probes probes = vfc_init_probes();
  probes = vfc_init_probes();

  bool ok;
  for (int i = 0; i < cycles_list.size(); i++) {
    ok = test_cycle(file, cycles_list[i], version, &probes);
    if (ok) {
      std::cout << "ok -- cycle " << std::to_string(i)
      << std::endl;
    }
    else {
      std::cerr << "failed -- cycle " << std::to_string(i)
      << std::endl;
    }
  }

  vfc_dump_probes(&probes);

  return ok;
}
