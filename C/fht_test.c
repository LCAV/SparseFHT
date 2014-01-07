
#include <stdio.h>
#include "fht.h"

int main(int argc, char **argv)
{
  int i;
  size_t len = 32;
  double X[32];

  for (i = 0 ; i < len/4 ; i++)
    X[i] = 1;
  for (i = len/4 ; i < len ; i++)
    X[i] = 0;

  /* allocate the memory */
  printf("The original vector: ");
  for (i = 0 ; i < len ; i++)
    printf("%f ", X[i]);
  printf("\n");

  fht(X, len);

  printf("The transformed vector: ");
  for (i = 0 ; i < len ; i++)
    printf("%f ", X[i]);
  printf("\n");

  return 0;
}
