#include "SMWB.hpp"
#include "Helpers.hpp"
#include "SM_Standard.hpp"
#include "Woodbury.hpp"

#define DEBUG1
// #define DEBUG2

// Sherman-Morrison-Woodbury kernel WB2s
// WB2, SM2 mixing scheme
void WB2s(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index, const double breakdown) {
#ifdef DEBUG2
  std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates
            << " updates" << std::endl;
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  unsigned int n_of_2blocks = N_updates / 2;
  unsigned int remainder = N_updates % 2;
  unsigned int length_2block = 2 * Dim;

  // Apply first 2*n_of_2blocks updates in n_of_2blocks blocks of 2 updates with
  // Woodbury 2x2 kernel
  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;
  if (n_of_2blocks > 0) {
    for (unsigned int i = 0; i < n_of_2blocks; i++) {
      double *Updates_2block = &Updates[i * length_2block];
      unsigned int *Updates_index_2block = &Updates_index[i * 2];
      bool ok;
      ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block, breakdown);
      if (!ok) { // Send the entire block to SM2
#ifdef DEBUG1
        std::cerr << "Woodbury 2x2 block failed! Sending to SM w/ US"
                  << std::endl;
#endif
#ifdef DEBUG2
        showMatrix2(Updates_2block, 2, Dim, "Updates_2block");
        showMatrix2(Updates_index_2block, 1, 2, "Updates_index_2block");
#endif
        unsigned int l = 0;
        SM2star(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block,
                later_updates + (Dim * later), later_index + later, &l, breakdown);
        later = later + l;
      }
    }
  }

  if (remainder != 0) { // Apply last remaining update with SM2
    double *Updates_1block = &Updates[n_of_2blocks * length_2block];
    unsigned int *Updates_index_1block = &Updates_index[2 * n_of_2blocks];
    unsigned int l = 0;
    SM2star(Slater_inv, Dim, remainder, Updates_1block, Updates_index_1block,
            later_updates + (Dim * later), later_index + later, &l, breakdown);
    later = later + l;
  }

  if (later > 0) {
    SM2(Slater_inv, Dim, later, later_updates, later_index, breakdown);
  }
}

// Sherman-Morrison-Woodbury kernel WB3s
// WB3, SM2 mixing scheme
void WB3s(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index, const double breakdown) {
#ifdef DEBUG2
  std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates
            << " updates" << std::endl;
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with
  // Woodbury 3x3 kernel
  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double *Updates_3block = &Updates[i * length_3block];
      unsigned int *Updates_index_3block = &Updates_index[i * 3];
      bool ok;
      ok = WB3(Slater_inv, Dim, Updates_3block, Updates_index_3block, breakdown);
      if (!ok) { // Send the entire block to SM2
#ifdef DEBUG1      
        std::cerr << "Woodbury 3x3 block failed! Sending to SM w/ US"
                  << std::endl;
#endif                  
#ifdef DEBUG2
        showMatrix2(Updates_3block, 3, Dim, "Updates_3block");
        showMatrix2(Updates_index_3block, 1, 3, "Updates_index_3block");
#endif
        unsigned int l = 0;
        SM2star(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block,
                later_updates + (Dim * later), later_index + later, &l, breakdown);
        later = later + l;
      }
    }
  }

  if (remainder != 0) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    unsigned int l = 0;
    SM2star(Slater_inv, Dim, remainder, Updates_2block, Updates_index_2block,
            later_updates + (Dim * later), later_index + later, &l, breakdown);
    later = later + l;
  }

  if (later > 0) {
    SM2(Slater_inv, Dim, later, later_updates, later_index, breakdown);
  }
}

// Sherman-Morrison-Woodbury kernel WB32s
// WB3, WB2, SM2 mixing scheme
void WB32s(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index, const double breakdown) {
#ifdef DEBUG2
  std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates
            << " updates" << std::endl;
  showMatrix2(Updates_index, 1, N_updates, "Updates_index");
  showMatrix2(Updates, N_updates, Dim, "Updates");
#endif

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with
  // Woodbury 3x3 kernel
  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double *Updates_3block = &Updates[i * length_3block];
      unsigned int *Updates_index_3block = &Updates_index[i * 3];
      bool ok;
      ok = WB3(Slater_inv, Dim, Updates_3block, Updates_index_3block, breakdown);
      if (!ok) { // Send the entire block to SM2
#ifdef DEBUG1      
        std::cerr << "Woodbury 3x3 block failed! Sending to SM w/ US"
                  << std::endl;
#endif                  
#ifdef DEBUG2
        showMatrix2(Updates_3block, 3, Dim, "Updates_3block");
        showMatrix2(Updates_index_3block, 1, 3, "Updates_index_3block");
#endif
        unsigned int l = 0;
        SM2star(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block,
                later_updates + (Dim * later), later_index + later, &l, breakdown);
        later = later + l;
      }
    }
  }

  if (remainder == 2) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    bool ok;
    ok = WB2(Slater_inv, Dim, Updates_2block, Updates_index_2block, breakdown);
    if (!ok) { // Send the entire block to SM2
#ifdef DEBUG1
      std::cerr << "Woodbury 2x2 block failed! Sending to SM w/ US"
                << std::endl;
#endif
      unsigned int l = 0;
      SM2star(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block,
              later_updates + (Dim * later), later_index + later, &l, breakdown);
      later = later + l;
    }
  } else if (remainder == 1) { // Apply last remaining update with SM2
    double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    unsigned int l = 0;
    SM2star(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block,
            later_updates + (Dim * later), later_index + later, &l, breakdown);
    later = later + l;
  }

  if (later > 0) {
    SM2(Slater_inv, Dim, later, later_updates, later_index, breakdown);
  }
}


extern "C" {
void WB2s_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
             double **linUpdates, unsigned int **Updates_index, const double breakdown) {
  WB2s(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index, breakdown);
}
void WB3s_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
             double **linUpdates, unsigned int **Updates_index, const double breakdown) {
  WB3s(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index, breakdown);
}
void WB32s_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
             double **linUpdates, unsigned int **Updates_index, const double breakdown) {
  WB32s(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index, breakdown);
}
}
