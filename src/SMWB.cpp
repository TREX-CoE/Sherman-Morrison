#include "SMWB.hpp"
#include "SM_Standard.hpp"
#include "Woodbury.hpp"
#include "Helpers.hpp"

// Sherman-Morrison-Woodbury kernel 1
// WB2, WB3, SM2 mixing scheme
void SMWB1(double *Slater_inv, unsigned int Dim, unsigned int N_updates, double *Updates, unsigned int *Updates_index) {
#ifdef DEBUG2
  std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates << " updates" << std::endl;
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with Woodbury 3x3 kernel
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double *Updates_3block = &Updates[i * length_3block];
      unsigned int *Updates_index_3block = &Updates_index[i * 3];
      bool ok;
      ok = WB3(Slater_inv, Dim, Updates_3block, Updates_index_3block);
      if (!ok) { // Send the entire block to SM2
#ifdef DEBUG2      
        std::cerr << "Woodbury 3x3 kernel failed! Sending block to SM2" << std::endl;
        showMatrix2(Updates_3block, 3, Dim, "Updates_3block");
        showMatrix2(Updates_index_3block, 1, 3, "Updates_index_3block");
#endif
        SM2(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block);
      }
    }
  }

  if (remainder == 2) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    bool ok;
    ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
    if (!ok) { // Send the entire block to SM2
#ifdef DEBUG2    
      std::cerr << "Woodbury 2x2 kernel failed! Sending block to SM2" << std::endl;
#endif      
      SM2(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block);
    }
  } else if (remainder == 1) { // Apply last remaining update with SM2
    double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    SM2(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block);
  } else { // remainder == 0
    // Nothing left to do.
  }
}

// Sherman-Morrison-Woodbury kernel 2
// WB2, WB3, SM3 mixing scheme
void SMWB2(double *Slater_inv, unsigned int Dim, unsigned int N_updates, double *Updates, unsigned int *Updates_index) {
#ifdef DEBUG2
  std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates << " updates" << std::endl;
#endif

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

#ifdef DEBUG2
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with Woodbury 3x3 kernel
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double *Updates_3block = &Updates[i * length_3block];
      unsigned int *Updates_index_3block = &Updates_index[i * 3];
      bool ok;
      ok = WB3(Slater_inv, Dim, Updates_3block, Updates_index_3block);
      if (!ok) { // Send the entire block to SM3
#ifdef DEBUG2      
        std::cerr << "Woodbury 3x3 kernel failed! Sending block to SM3" << std::endl;
        showMatrix2(Updates_3block, 3, Dim, "Updates_3block");
        showMatrix2(Updates_index_3block, 1, 3, "Updates_index_3block");
#endif
        SM3(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block);
      }
    }
  }

  if (remainder == 2) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    bool ok;
    ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
    if (!ok) { // Send the entire block to SM3
      std::cerr << "Woodbury 2x2 kernel failed! Sending block to SM3" << std::endl;
      SM3(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block);
    }
  } else if (remainder == 1) { // Apply last remaining update with SM3
    double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    SM3(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block);
  } else { // remainder == 0
    // Nothing left to do.
  }
}

// Sherman-Morrison-Woodbury kernel 3
// WB2, WB3, SM4 mixing scheme
void SMWB3(double *Slater_inv, unsigned int Dim, unsigned int N_updates, double *Updates, unsigned int *Updates_index) {
           std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates << " updates" << std::endl;

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

#ifdef DEBUG2
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with Woodbury 3x3 kernel
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double *Updates_3block = &Updates[i * length_3block];
      unsigned int *Updates_index_3block = &Updates_index[i * 3];
      bool ok;
      ok = WB3(Slater_inv, Dim, Updates_3block, Updates_index_3block);
      if (!ok) { // Send the entire block to SM4
        std::cerr << "Woodbury 3x3 kernel failed! Sending block to SM4" << std::endl;
#ifdef DEBUG2        
        showMatrix2(Updates_3block, 3, Dim, "Updates_3block");
        showMatrix2(Updates_index_3block, 1, 3, "Updates_index_3block");
#endif
        SM4(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block);
      }
    }
  }

  if (remainder == 2) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    bool ok;
    ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
    if (!ok) { // Send the entire block to SM4
      std::cerr << "Woodbury 2x2 kernel failed! Sending block to SM4" << std::endl;
      SM4(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block);
    }
  } else if (remainder == 1) { // Apply last remaining update with SM4
    double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    SM4(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block);
  } else { // remainder == 0
    // Nothing left to do.
  }
}

// Sherman-Morrison-Woodbury kernel 4
// WB2, SM2 mixing scheme
void SMWB4(double *Slater_inv, unsigned int Dim, unsigned int N_updates, double *Updates, unsigned int *Updates_index) {
#ifdef DEBUG2
  std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates << " updates" << std::endl;
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  unsigned int n_of_2blocks = N_updates / 2;
  unsigned int remainder = N_updates % 2;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

  // Apply first 2*n_of_2blocks updates in n_of_2blocks blocks of 2 updates with Woodbury 2x2 kernel
  if (n_of_2blocks > 0) {
    for (unsigned int i = 0; i < n_of_2blocks; i++) {
      double *Updates_2block = &Updates[i * length_2block];
      unsigned int *Updates_index_2block = &Updates_index[i * 2];
      bool ok;
      ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
      if (!ok) { // Send the entire block to SM2
        std::cerr << "Woodbury 2x2 kernel failed! Sending block to SM2" << std::endl;
#ifdef DEBUG2
        showMatrix2(Updates_2block, 2, Dim, "Updates_2block");
        showMatrix2(Updates_index_2block, 1, 2, "Updates_index_2block");
#endif
        SM2(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block);
      }
    }
  }

  if (remainder == 1) { // Apply last remaining update with SM4
    double *Updates_1block = &Updates[n_of_2blocks * length_2block];
    unsigned int *Updates_index_1block = &Updates_index[2 * n_of_2blocks];
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
void SMWB2_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
           double **linUpdates, unsigned int **Updates_index) {
  SMWB2(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index);
}
void SMWB3_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
           double **linUpdates, unsigned int **Updates_index) {
  SMWB3(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index);
}
void SMWB4_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
           double **linUpdates, unsigned int **Updates_index) {
  SMWB4(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index);
}
}
