// Woodbury 2x2 kernel
bool WB2(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index, const double breakdown);

// Woodbury 3x3 kernel
bool WB3(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index, const double breakdown);
