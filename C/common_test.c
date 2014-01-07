
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"

#define N 10
#define MAX 1024

int main(int argc, char **argv)
{

  int i;
  int array[N];

  srandom(time(NULL));

  printf("Start: ");
  for (i = 0 ; i < N ; i++)
  {
    array[i] = random() % MAX;
    printf("%d ", array[i]);
  }
  printf("\n");

  quicksort(array, 0, N-1);

  printf("Result: ");
  for (i = 0 ; i < N ; i++)
    printf("%d ", array[i]);
  printf("\n");

  return 0;
}
