
#include "fht.h"

/* a simple recursive in-place fast hadamard transform */
/* x is the input vector, and n is the size, a power of two */
void fht(double *x, size_t n)
{
  /* stopping case */
  if (n == 1)
    return;

  /* recurse */
  fht(x, n/2);
  fht(x+n/2, n/2);

  /* butterflies */
  int i;
  for (i = 0 ; i < n/2 ; i++)
  {
    double tmp = x[i];
    x[i] = tmp + x[n/2+i];
    x[n/2+i] = tmp - x[n/2+i];
  }

  return;
}
