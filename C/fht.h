#ifndef __FHT_H__
#define __FHT_H__

#include <stdlib.h>

/* a simple recursive in-place fast hadamard transform */
/* x is the input vector, and n is the size, a power of two */
void fht(double *x, size_t n);

#endif /* __FHT_H__ */
