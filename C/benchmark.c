
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "fht.h"
#include "sfht.h"

#define N 32768
#define K 4096
#define ALPHA 2
#define C 2
#define L 2
#define MAX 512
#define WARMUP 10
#define LOOP 20

int main(int argc, char **argv)
{
  int i, q;
  double x[N];
  int index[N];
  double y[K];
  int s[K];
  time_t now;

  unsigned seed = 0;
  seed = time(NULL);

  srandom(seed);

  printf("Benchmark -- start:\n");

  for (i = 0 ; i < N ; i++)
  {
    index[i] = i;
    x[i] = 0;
  }

  for (i = 0 ; i < K ; i++)
  {
    int p = (random() % N-i) + i;
    int t = index[p];
    index[p] = index[i];
    index[i] = t;
  }

  printf("SparseFHT... ");
  for (i = 0 ; i < K ; i++)
    x[index[i]] = (random() % MAX) * (1 - 2*random()%2);

  for (q = 0 ; q < WARMUP ; q++)
    sfht_deterministic(x, N, y, s, K, K, C, L, NULL, NULL);
  //sfht(x, N, y, s, K, ALPHA, C, L);
  now = clock();
  for (q = 0 ; q < LOOP ; q++)
    sfht_deterministic(x, N, y, s, K, K, C, L, NULL, NULL);
    //sfht(x, N, y, s, K, ALPHA, C, L);
  now = clock() - now;
  printf("Avg runtime: %u \n", (unsigned)(now/LOOP));
  fflush(stdout);

  printf("FHT...       ");
  for (q = 0 ; q < WARMUP ; q++)
    fht(x, N);
  now = clock();
  for (q = 0 ; q < LOOP ; q++)
    fht(x, N);
  now = clock() - now;
  printf("Avg runtime: %u \n", (unsigned)(now/LOOP));

  return 0;
}
