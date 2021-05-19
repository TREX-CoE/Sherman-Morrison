#include "SM_Helpers.hpp"

// Set common break-down threshold
double threshold() {
  const double threshold = THRESHOLD;
#ifdef DEBUG
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
#ifdef DEBUG
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
