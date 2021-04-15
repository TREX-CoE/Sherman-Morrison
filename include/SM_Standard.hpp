// Naïve Sherman Morrison
void SM1(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index);

// Sherman Morrison, with J. Slagel splitting
// http://hdl.handle.net/10919/52966
void SM2(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index);

// Sherman Morrison, leaving zero denominators for later
void SM3(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index);
