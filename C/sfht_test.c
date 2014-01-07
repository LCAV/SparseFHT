
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "fht.h"
#include "sfht.h"
#include "common.h"

#define N (1 << 19)
#define K (1 << 12)
#define MAX 25
#define LOOP 100
#define ALPHA 1
#define C 3
#define L 10

int main(int argc, char **argv)
{
  int i, q;
  double x[N];
  gf2_vector_t index[N];
  double y[K];
  int s[K];
  double avg_error = 0;
  int loops, unsat;
  double n_error = 0;

  srandom(time(NULL));

  printf("Test Sparse Hadamard :\n");

  for (i = 0 ; i < N ; i++)
    index[i] = i;

  for (q = 0 ; q < LOOP ; q++)
  {
    /* random permutation */
    for (i = 0 ; i < N-1 ; i++)
    {
      int p = random() % (N-i) + i;
      double t = index[p];
      index[p] = index[i];
      index[i] = t;
    }

    /* create sparse vector */
    for (i = 0 ; i < N ; i++)
      x[i] = 0;
    for (i = 0 ; i < K ; i++)
    {
      x[index[i]] = random() % (MAX-1) + 1;
      //s[i] = index[i];
      //y[i] = x[index[i]];
    }

    /* prepare TD signal */
    fht(x, N);
    for (i = 0 ; i < N ; i++)
      x[i] /= sqrt(N);

    /* sparse transform */
    srandom(100);
    sfht_random(x, N, y, s, K, ALPHA*K, C, L, &loops, &unsat);
    srandom(time(NULL));
    
    if (unsat > 0)
      n_error++;

    /* recover FD signal */
    fht(x, N);
    for (i = 0 ; i < N ; i++)
      x[i] /= sqrt(N);

    sparse_quicksort(y, (gf2_vector_t *)s, 0, K-1);

    int u = 0;
    int v = 0;
    double error = 0.0;
    quicksort((int *)index, 0, K-1);
    while (u < K && v < K)
    {
      while (u < K && (v >= K || s[u] < index[v]))
        error += fabs(y[u++]);
      while (v < K && (u >= K || index[v] < s[u]))
        error += fabs(x[index[v++]]);
      u++;
      v++;
    }
    printf("Error: %.3f (loops=%d, unsat=%d)\n", error, loops, unsat);
    avg_error += error;
      
  }
  printf("Average Error: %.3f\n", avg_error/LOOP);
  printf("Failure rate: %.3f\n", n_error/(float)LOOP);

#if 0 
  if (LOOP == 1)
  {
    for (i = 0 ; i < K ; i++)
    {
      printf("O(%d,%.2f) S(%d,%.2f) ", index[i], x[index[i]], s[i], y[i]);
      if (fabs(x[index[i]] - y[i]) > 1e-5)
        printf(" <---- check this");
      printf("\n");
    }
  }
#endif

  return 0;
}
