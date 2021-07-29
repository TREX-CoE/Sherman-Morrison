// Sherman-Morrison-Woodbury kernel WB2s
// WB2, SM2 mixing scheme
void WB2s(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index, const double breakdown);

// Sherman-Morrison-Woodbury kernel WB3s
// WB3, SM2 mixing scheme
void WB3s(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index, const double breakdown);

// Sherman-Morrison-Woodbury kernel WB32s
// WB3, WB2, SM2 mixing scheme
void WB32s(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index, const double breakdown);
