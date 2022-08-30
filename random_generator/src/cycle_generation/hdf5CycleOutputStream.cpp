#include "cycle_generation/hdf5CycleOutputStream.hpp"

namespace fs = std::filesystem;

namespace randomgen {


  namespace {

    void writeDeterminant(const Cycle &cycle,
                          hid_t cycle_gid) {// Write the determinant of the matrix
      hid_t det_space = H5Screate_simple(1, std::vector<hsize_t>{1}.data(), nullptr);
      hid_t determinant_loc =
              H5Dcreate1(cycle_gid, "determinant", H5T_NATIVE_DOUBLE, det_space, H5P_DEFAULT);
      if (determinant_loc < 0) { throw std::runtime_error("H5Dcreate1"); }

      double determinant = cycle.getSlaterMatrix().getLastKnownDeterminant();

      H5Dwrite(determinant_loc, H5T_NATIVE_DOUBLE, det_space, H5S_ALL, H5S_ALL, &determinant);
      H5Dclose(determinant_loc);
      H5Sclose(det_space);
    }

    void writeMatrix(const Matrix &mat, hid_t group_id, hid_t slater_space,
                     const std::string &dataset_name) {
      hid_t slater_matrix_loc = H5Dcreate1(group_id, dataset_name.c_str(), H5T_NATIVE_DOUBLE,
                                           slater_space, H5P_DEFAULT);
      if (slater_matrix_loc < 0) { throw std::runtime_error("H5Dcreate1"); }

      H5Dwrite(slater_matrix_loc, H5T_NATIVE_DOUBLE, slater_space, H5S_ALL, H5S_ALL, mat.data());
      H5Dclose(slater_matrix_loc);
    }


    void writeMatrices(const Cycle &cycle, hid_t group_gid) {
      unsigned int y_size = cycle.getSlaterInverseTransposed().y,
                   x_size = cycle.getSlaterInverseTransposed().x;

      // Write both the matrix and its inverse
      auto slater_space = H5Screate_simple(2, std::vector<hsize_t>{x_size, y_size}.data(), nullptr);
      writeMatrix(cycle.getSlaterMatrix(), group_gid, slater_space, "slater_matrix");
      writeMatrix(cycle.getSlaterInverseTransposed(), group_gid, slater_space, "slater_inverse_t");
      H5Sclose(slater_space);

      writeDeterminant(cycle, group_gid);
    }

    void writeUpdatesMetadata(const Cycle &cycle, hid_t cycle_gid) {// Write the number of updates
      hid_t updates_count_space = H5Screate_simple(1, std::vector<hsize_t>{1}.data(), nullptr);
      hid_t updates_count_loc = H5Dcreate1(cycle_gid, "nupdates", H5T_NATIVE_UINT32,
                                           updates_count_space, H5P_DEFAULT);
      if (updates_count_loc < 0) { throw std::runtime_error("H5Dcreate1"); }

      auto update_count = cycle.getUpdateMatrix().getUpdateCount();
      H5Dwrite(updates_count_loc, H5T_NATIVE_UINT32, updates_count_space, H5S_ALL, H5S_ALL,
               &update_count);
      H5Dclose(updates_count_loc);
      H5Sclose(updates_count_space);


      // Write the columns of the updates matrix
      hid_t updates_index = H5Screate_simple(
              1, std::vector<hsize_t>{cycle.getUpdateMatrix().getUpdateCount()}.data(), nullptr);
      hid_t updates_index_loc = H5Dcreate1(cycle_gid, "col_update_index", H5T_NATIVE_UINT32,
                                           updates_index, H5P_DEFAULT);
      if (updates_index_loc < 0) { throw std::runtime_error("H5Dcreate1"); }
      H5Dwrite(updates_index_loc, H5T_NATIVE_UINT32, updates_index, H5S_ALL, H5S_ALL,
               cycle.getUpdateMatrix().getUpdateIndices());
      H5Dclose(updates_index_loc);
      H5Sclose(updates_index);
    }

    void writeUpdateMatrix(const Cycle &cycle, hid_t cycle_gid) {// Write the updates matrix
      unsigned int y_size = cycle.getSlaterInverseTransposed().y;

      hid_t update_space = H5Screate_simple(
              2, std::vector<hsize_t>{cycle.getUpdateMatrix().getUpdateCount(), y_size}.data(),
              nullptr);
      hid_t update_loc =
              H5Dcreate1(cycle_gid, "updates", H5T_NATIVE_DOUBLE, update_space, H5P_DEFAULT);
      if (update_loc < 0) { throw std::runtime_error("H5Dcreate1"); }
      H5Dwrite(update_loc, H5T_NATIVE_DOUBLE, update_space, H5S_ALL, H5S_ALL,
               cycle.getUpdateMatrix().getRawUpdates());
      H5Dclose(update_loc);
      H5Sclose(update_space);
    }

    void writeUpdates(const Cycle &cycle, hid_t cycle_id) {
      writeUpdatesMetadata(cycle, cycle_id);
      writeUpdateMatrix(cycle, cycle_id);
    }

  }// namespace

  hdf5CycleOutputStream::hdf5CycleOutputStream(const fs::path &dataset_path) {
    // Create a new file
    file_id = H5Fcreate(dataset_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) { throw std::runtime_error("Could not create dataset file"); }
  }

  hdf5CycleOutputStream::~hdf5CycleOutputStream() {
    if (file_id >= 0) { H5Fclose(file_id); }
    file_id = -1;
  }

  hid_t hdf5CycleOutputStream::createGroup(const std::string &path) {
    // Needed to create the parent directories automatically
    auto lcpl = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(lcpl, 1);
    hid_t group_id = H5Gcreate2(file_id, path.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(lcpl);
    return group_id;
  }

  std::string hdf5CycleOutputStream::getPathFor(const Cycle &cycle) {
    std::string res = "/";
    res += std::to_string(cycle.getUpdateMatrix().getUpdateCount()) + "/";
    res += std::to_string(cycle.getUpdateMatrix().getSplitCount()) + "/";
    res += "cycle_" + std::to_string(last_unique_cycle_id);
    return res;
  }

  hdf5CycleOutputStream &operator<<(hdf5CycleOutputStream &ds, const Cycle &cycle) {
    // Build the path to the new dataset inside the file
    std::string cycle_path = ds.getPathFor(cycle);
    hid_t group_id = ds.createGroup(cycle_path);
    ds.last_unique_cycle_id++;

    writeMatrices(cycle, group_id);
    writeUpdates(cycle, group_id);

    H5Gclose(group_id);
    return ds;
  }
}// namespace randomgen
