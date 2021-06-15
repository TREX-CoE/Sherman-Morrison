#include "SMWB.hpp"
#include "SM_Standard.hpp"
#include "Woodbury.hpp"
#include "Helpers.hpp"

// Sherman-Morrison-Woodbury kernel 1
// WB2, WB3, SM2 mixing scheme 1
void SMWB1(double *Slater_inv, unsigned int Dim, unsigned int N_updates, double *Updates, unsigned int *Updates_index) {
           std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates << " updates" << std::endl;

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

  // std::cerr << "Number of blocks: " << n_of_3blocks << ". Remainder: " << remainder << "." << std::endl;

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with Woodbury 3x3 kernel
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double Updates_3block[length_3block];
      unsigned int Updates_index_3block[3];
      Updates_index_3block[0] = Updates_index[3 * i + 0];
      Updates_index_3block[1] = Updates_index[3 * i + 1];
      Updates_index_3block[2] = Updates_index[3 * i + 2];
      for (unsigned int j = 0; j < length_3block; j++) {
        Updates_3block[j] = Updates[i * length_3block + j];
      }
      bool ok;
      ok = WB3(Slater_inv, Dim, Updates_3block, Updates_index_3block);
      if (!ok) { // Send the entire block to SM2
        std::cerr << "Woodbury 3x3 kernel failed! Sending block to SM2" << std::endl;
        SM2(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block);
      }
    }
  }

  if (remainder == 2) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double Updates_2block[length_2block];
    unsigned int Updates_index_2block[2];
    Updates_index_2block[0] = Updates_index[3 * n_of_3blocks + 0];
    Updates_index_2block[1] = Updates_index[3 * n_of_3blocks + 1];
    for (unsigned int i = 0; i < length_2block; i++) {
      Updates_2block[i] = Updates[n_of_3blocks * length_3block + i];
    }
    bool ok;
    ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
    if (!ok) { // Send the entire block to SM2
      std::cerr << "Woodbury 2x2 kernel failed! Sending block to SM2" << std::endl;
      SM2(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block);
    }
  } else if (remainder == 1) { // Apply last remaining update with SM2
    double Updates_1block[length_1block];
    unsigned int Updates_index_1block[1];
    Updates_index_1block[0] = Updates_index[3 * n_of_3blocks + 0];
    for (unsigned int i = 0; i < length_1block; i++) {
      Updates_1block[i] = Updates[n_of_3blocks * length_3block + i];
    }
    SM2(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block);
  } else { // remainder == 0
    // Nothing left to do.
  }
}

extern "C" {
void SMWB1_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
           double **linUpdates, unsigned int **Updates_index) {
  SMWB1(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index);
}
}
