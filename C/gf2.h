#ifndef __GF2_H__
#define __GF2_H__
/*
 * Implement a few matrix operations in GF(2)
 *
 * Every row of a GF(2) matrix is represented by a 32 bit integer
 * thus limiting the size of the matrices to 32x32.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* maximum size allowed by underlying 32 bit words */
/* XXX For some reasong there is a bug for 32, set to 31 for now XXX */
#define GF2_MAX_SIZE 31

/* just a wrapper for binary vectors of at most 32 bits */
typedef uint32_t gf2_vector_t;

/* 
 * A square GF(2) matrix of size NxN
 * U is a list of integer such that
 * The N right-most bits of U[r] are
 * the elements of the r^th row of U.
 */
typedef struct {
  gf2_vector_t *rows;
  size_t size;
} gf2_matrix_t;

/* 
 * Create a random non-singular matrix U in GF(2) and its inverse V.
 */ 
gf2_matrix_t *gf2_rand_mat(gf2_matrix_t *U, gf2_matrix_t *V);

/*
 * Generate random non-singular matrix efficiently
 * Based on Dana Randall's paper
 */
void gf2_rand_nonsing(gf2_matrix_t *U);

/* allocate and free GF(2) matrices */
gf2_matrix_t *gf2_alloc_mat(size_t N);
void gf2_free_mat(gf2_matrix_t *U);

/* computes the inverse of U in GF(2) and stores it in V */
gf2_matrix_t *gf2_inverse(gf2_matrix_t *U, gf2_matrix_t *V);

/* Transposes U */
void gf2_transpose(gf2_matrix_t *U);

/* 
 * multiplies a vector V in GF(2) stored in an integer by matrix U infers the
 * length of V from the size of U 
 */
gf2_vector_t gf2_mat_vec(gf2_matrix_t *U, gf2_vector_t V);

/*
 * compute product of matrix and a number of vectors, in-place
 */
void gf2_mat_vec_array(gf2_matrix_t *U, gf2_vector_t *V, size_t len);

/*
 * compute the number of ones in the vector V of length N
 */
#if 0
gf2_vector_t gf2_parity_32(gf2_vector_t V);
gf2_vector_t gf2_parity_16(gf2_vector_t V);
gf2_vector_t gf2_parity_8(gf2_vector_t V);
#endif
int gf2_parity(gf2_vector_t V);

/*
 * copy matrix U into V
 */
void gf2_copy(gf2_matrix_t *U, gf2_matrix_t *V);

/* simple printing routines */
void gf2_print_vec(gf2_vector_t V, size_t N);
void gf2_print_mat(gf2_matrix_t *U);

#endif /* __GF2_H__ */
