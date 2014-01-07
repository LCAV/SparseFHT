
#include <stdio.h>
#include <time.h>
#include "gf2.h"

int main(int argc, char **argv)
{
  int i;
  gf2_matrix_t *U, *V;
  size_t N = 5;
  int len = 1 << N;
  gf2_vector_t X[1<< N];

  for (i = 0 ; i < len ; i++)
    X[i] = i;

  srandom(time(NULL));

  /* allocate the memory */
  printf("Allocate the matrices.\n");
  U = gf2_alloc_mat(N);
  V = gf2_alloc_mat(N);

#if 0
  printf("Fill random matrix and compute inverse.\n");
  if (gf2_rand_mat(U, V) == NULL)
  {
    fprintf(stderr, "Error: could not allocate matrices U and V.\n");
    return 1;
  }
#endif
  gf2_rand_nonsing(U);
  if (gf2_inverse(U,V) == NULL)
    printf("The inverse is singular.\n");


  printf("The matrix: \n");
  gf2_print_mat(U);

  printf("The inverse: \n");
  gf2_print_mat(V);

  printf("The original vector: ");
  for (i = 0 ; i < len ; i++)
    printf("%d ", X[i]);
  printf("\n");

  gf2_mat_vec_array(U, X, len);

  printf("The permutation vector: ");
  for (i = 0 ; i < len ; i++)
    printf("%d ", X[i]);
  printf("\n");

  gf2_mat_vec_array(V, X, len);

  printf("And back: ");
  for (i = 0 ; i < len ; i++)
    printf("%d ", X[i]);
  printf("\n");

  printf("The inverse transpose: \n");
  gf2_transpose(V);
  gf2_print_mat(V);
  
  gf2_free_mat(U);
  gf2_free_mat(V);

  return 0;
}
