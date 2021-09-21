// Naïve Sherman Morrison
void SM1(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index, const double breakdown);

// Sherman Morrison, with J. Slagel splitting
// http://hdl.handle.net/10919/52966
void SM2(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index, const double breakdown);

void SM2star(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
             double *Updates, unsigned int *Updates_index,
             double *later_updates, unsigned int *later_index,
             unsigned int *later, const double breakdown);

// Sherman Morrison, leaving zero denominators for later
void SM3(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index, const double breakdown);

// Sherman Morrison (SM3+SM2), leaving zero denominators for later (SM3), and
// when none are left falling back on Splitting (SM2)
void SM4(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index, const double breakdown);
