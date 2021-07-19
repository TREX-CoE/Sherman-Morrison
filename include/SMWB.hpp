// Sherman-Morrison-Woodbury kernel 1
// WB2, WB3, SM2 mixing scheme
void SMWB1(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index);

// Sherman-Morrison-Woodbury kernel 2
// WB2, SM2 mixing scheme
void SMWB2(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index);
