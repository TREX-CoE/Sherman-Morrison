#include "SMWB.hpp"
#include "SM_Standard.hpp"
#include "Woodbury.hpp"
#include "Helpers.hpp"

// Sherman-Morrison-Woodbury kernel 1
// WB2, WB3, SM2 mixing scheme 1
void SMWB1(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index) {
           std::cerr << "Called Sherman-Morrison-Woodbury kernel 1 with " << N_updates << " updates" << std::endl;
           WB3(Slater_inv, Dim, Updates, Updates_index);
           WB2(Slater_inv, Dim, Updates, Updates_index);
           SM2(Slater_inv, Dim, N_updates, Updates, Updates_index);
         }

extern "C" {
void SMWB1_f(double **linSlater_inv, unsigned int *Dim, unsigned int *N_updates,
           double **linUpdates, unsigned int **Updates_index) {
  SMWB1(*linSlater_inv, *Dim, *N_updates, *linUpdates, *Updates_index);
}
}
