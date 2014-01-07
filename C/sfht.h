#ifndef __SFHT_H__
#define __SFHT_H__

#include <stdlib.h>
#include "gf2.h"

/* SFHT with random hash matrix */
int sfht_random(double *x, size_t N, double *y, int *s, size_t K, size_t B, int C, int L, int *loops_done, int *unsat);

/* SFHT with deterministic hash matrix */
int sfht_deterministic(double *x, size_t N, double *y, int *s, size_t K, size_t B, int C, int L, int *loops_done, int *unsat);

/* Underlying SFHT algorithm, hash matrix agnostic */
int sfht_base(double *x, size_t N, double *y, int *s, size_t K, 
    size_t B, int C, gf2_matrix_t **U, gf2_matrix_t **V, gf2_matrix_t **Vinv, int L, int *loops_done, int *unsat);

/* Sparse Hadamard hashing function */
void fast_hadamard_hash(double *x, size_t n, double *hash, size_t B, gf2_matrix_t *sigma, gf2_vector_t a);

/* 
 * Sparse Fast Hadamard Transform
 * using deterministic sampling pattern
 * and optimized bitwise operations to perform them
 */
int sfht_det_opt(double *x, size_t N, double *y, int *s, size_t K, size_t B, int C, int L, int *loops_done, int *unsat);

/* the deterministic sampling using bitwise operations */
uint32_t select_hash(uint32_t b, uint32_t n, uint32_t c, uint32_t C);
uint32_t hashing(uint32_t val, uint32_t a, uint32_t b, uint32_t n);
uint32_t inv_hashing(uint32_t val, uint32_t a, uint32_t b, uint32_t n);
uint32_t null_vector(uint32_t d, uint32_t a, uint32_t b, uint32_t n);

#endif /* __SFHT_H__ */
