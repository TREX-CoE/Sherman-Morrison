#include "Helpers.hpp"

// Set common break-down threshold
double threshold() {
  const double threshold = THRESHOLD;
#ifdef DEBUG2
  std::cerr << "Break-down threshold set to: " << threshold << std::endl;
#endif
  return threshold;
}

void Switch(unsigned int *p, unsigned int l, unsigned int lbar) {
  unsigned int tmp = p[l + 1];
  p[l + 1] = p[lbar];
  p[lbar] = tmp;
}

void selectLargestDenominator(unsigned int l, unsigned int N_updates,
                              unsigned int *Updates_index, unsigned int *p,
                              double ***ylk) {
  unsigned int lbar = l + 1, max = 0;
  unsigned int index = 0, component = 0;
  unsigned int tmp = 0;
  double breakdown = 0;
  for (unsigned int j = lbar; j < N_updates + 1; j++) {
    index = p[j];
    component = Updates_index[index - 1];
    breakdown = std::fabs(1 + ylk[l][index][component]);
#ifdef DEBUG2
    std::cout << "Inside selectLargestDenominator()" << std::endl;
    std::cout << "breakdown = fabs(1 + ylk[" << l << "][" << index << "]["
              << component << "]) = " << breakdown << std::endl;
    std::cout << std::endl;
#endif
    if (breakdown > max) {
      max = breakdown;
      lbar = j;
    }
  }
  Switch(p, l, lbar);
}

#ifdef MKL
// Inplace inverse n x n matrix A.
// returns:
//   ret = 0 on success
//   ret < 0 illegal argument value
//   ret > 0 singular matrix
lapack_int inverse(double *A, unsigned n) {
  int ipiv[n + 1];
  lapack_int ret;

  ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);

  if (ret != 0)
    return ret;

  ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, A, n, ipiv);
  return ret;
}
#endif
