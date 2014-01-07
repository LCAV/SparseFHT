#ifndef __COMMON_H__
#define __COMMON_H__

#include "gf2.h"

/* a simple quicksort */
void quicksort(int *array, int l, int r);

/* sparse array quicksort */
void swap2(double* value, gf2_vector_t *index, int k, int l);
void sparse_quicksort(double *value, gf2_vector_t *index, int l, int r);

#endif /* __COMMON_H__ */
