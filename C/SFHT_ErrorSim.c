
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "fht.h"
#include "sfht.h"
#include "common.h"

#define LOOP 10
#define MAX (1<<30)
#define L 10
#define SEED time(NULL)

int main(int argc, char **argv)
{
  size_t n = 22;
  size_t BLEN = 20;
  size_t b[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  size_t C[] = {7, 6, 5, 4, 3, 3, 3, 3, 3,  3,  3,  3,  3,  3,  4,  4,  5,  6,  8, 11};

  double error;
  double biterror;
  double n_loops;
  double success;

  size_t N = (1 << n);

  int i, q, k;
  double *x = NULL;
  gf2_vector_t *index = NULL;
  double *y = NULL;
  int *s = NULL;
  int loops, unsat;

  if (argc != 2)
    printf("Usage: %s <file>\n", argv[0]);

  FILE *outfile = fopen(argv[1], "w");
  if (outfile == NULL)
  {
    fprintf(stderr, "File %s can't open output file.\n", argv[1]);
    exit(1);
  }

  srandom(SEED);

  printf("Test Sparse Hadamard :\n");

  x = (double *)malloc(N*sizeof(double));
  index = (gf2_vector_t *)malloc(N*sizeof(gf2_vector_t));
  for (i = 0 ; i < N ; i++)
    index[i] = i;

  for (k = 0 ; k < BLEN ; k++)
  {
    /* init trackers */
    n_loops = 0;
    success = 0;
    error = 0;
    biterror = 0;

    size_t K = (1 << b[k]);

    /* allocate arrays */
    y = (double *)malloc(K*sizeof(double));
    s = (int *)malloc(K*sizeof(int));
    if (x == NULL || index == NULL || y == NULL || s == NULL)
    {
      fprintf(stderr, "Could not allocate memory.\n");
      exit(1);
    }


    for (q = 0 ; q < LOOP ; q++)
    {
      /* random permutation */
      for (i = 0 ; i < K ; i++)
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
        x[index[i]] = ((random()%MAX)/(float)MAX)-0.5;
        s[i] = index[i];
        y[i] = x[index[i]];
      }

      /* prepare TD signal */
      fht(x, N);
      for (i = 0 ; i < N ; i++)
        x[i] /= sqrt(N);

      /* sparse transform */
      //int ncoeff = sfht_det_opt(x, (size_t)N, y, s, K, (1 << b[k]), C[k], L, &loops, &unsat);
      int ncoeff = sfht_deterministic(x, (size_t)N, y, s, K, (1 << b[k]), C[k], L, &loops, &unsat);
    
      
      n_loops += loops;
      //if (unsat <= 0)
       //uccess++;

      /* recover FD signal */
      fht(x, N);
      for (i = 0 ; i < N ; i++)
        x[i] /= sqrt(N);

      int u = 0;
      int local_error = 0;

      for (u = 0 ; u < ncoeff ; u++)
        x[s[u]] -= y[u];

      for (u = 0 ; u < N ; u++)
        if (fabs(x[u]) > 1e-5)
        {
          local_error += 1;
          error += fabs(x[u]);
          biterror += 1;
        }

      if (local_error == 0)
        success += 1;
      //printf("%f %f %f %f %d %d\n", error, biterror, n_loops, success, loops, unsat);

    }
    printf("b=%d alpha=%.3f Error: %f BitError: %f (loops=%f, success=%f)\n", (int)b[k], b[k]/(float)n, error/LOOP, biterror/LOOP, n_loops/LOOP, success/LOOP);
    fprintf(outfile, "%f %f %f %f %d %d %d %d\n", error, biterror, n_loops, success, (int)n, (int)b[k], (int)C[k], (int)LOOP);
    //fprintf(outfile, "hello %d\n", k);
        
    free(y);
    free(s);

  }

  fclose(outfile);

  /* free mem */
  free(x);
  free(index);

  return 0;
}
