
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "fht.h"
#include "sfht.h"
#include "common.h"

/*
#define N 524288
#define K 16384
#define MAX 25
#define LOOP 100
#define C 3
#define L 10
#define SEED time(NULL)
*/

#define N (1 << 22)
#define K (1 << 18)
#define MAX 1000000l
#define LOOP 100
#define C 11
#define L 10
#define SEED time(NULL)

/*
#define N 524288
#define K 8
#define MAX 25
#define LOOP 100
#define C 3
#define L 10
#define SEED time(NULL)
*/

int main(int argc, char **argv)
{
  int i, q;
  double *x;
  gf2_vector_t *index;
  double *y;
  int *s;
  double avg_error=0;
  int loops, unsat;
  int n_error = 0;

  /* allocate arrays */
  x = (double *)malloc(N*sizeof(double));
  index = (gf2_vector_t *)malloc(N*sizeof(gf2_vector_t));
  y = (double *)malloc(K*sizeof(double));
  s = (int *)malloc(K*sizeof(int));
  if (x == NULL || index == NULL || y == NULL || s == NULL)
  {
    fprintf(stderr, "Could not allocate memory.\n");
    exit(1);
  }

  srandom(SEED);

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
      x[index[i]] = ((random()%MAX)/(float)MAX);
      s[i] = index[i];
      y[i] = x[index[i]];
    }

    /* prepare TD signal */
    fht(x, N);
    for (i = 0 ; i < N ; i++)
      x[i] /= sqrt(N);

    /* sparse transform */
    int ncoeff = sfht_det_opt(x, (size_t)N, y, s, K, K, C, L, &loops, &unsat);
    //int ncoeff = sfht_deterministic(x, (size_t)N, y, s, K, K, C, L, &loops, &unsat);
    
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
    int biterror = 0;
    quicksort((int *)index, 0, K-1);
    while (u < K && v < K)
    {
      while (u < K && (v >= K || s[u] < index[v]))
      {
        error += fabs(y[u++]);
        biterror++;
      }
      while (v < K && (u >= K || index[v] < s[u]))
      {
        error += fabs(x[index[v++]]);
        biterror++;
      }
      u++;
      v++;
    }
    printf("Error: %.3f BitError: %d (ncoeff=%d, loops=%d, unsat=%d)\n", error, biterror, ncoeff, loops, unsat);
    avg_error += error;
      
  }
  printf("Average Error: %.3f\n", avg_error/LOOP);
  printf("Error rate: %.3f\n", n_error/(float)LOOP);

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

  /* free mem */
  free(x);
  free(index);
  free(y);
  free(s);

  return 0;
}
