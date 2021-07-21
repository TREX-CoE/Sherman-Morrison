#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#ifndef THRESHOLD
#define THRESHOLD 1e-3
#endif

int main() {
  qmckl_context context;
  context = qmckl_context_create();

  qmckl_exit_code rc;

const uint64_t Dim = 2;
const uint64_t N_updates = 2;
const uint64_t Updates_index[2] = {0, 0};
const double Updates[4] = {0.0, 0.0, 0.0, 0.0};
double Slater_inv[4] = {0.0, 0.0, 0.0, 0.0};

rc = qmckl_sherman_morrison_c(context, Dim, N_updates, Updates, Updates_index, Slater_inv);
assert(rc == QMCKL_SUCCESS);

assert (qmckl_context_destroy(context) == QMCKL_SUCCESS);
  return 0;
}
