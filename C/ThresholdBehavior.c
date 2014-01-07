
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "fht.h"
#include "sfht.h"
#include "common.h"

#define ALEN 20
#define MAX 1000000l
#define LOOP 100
#define L 10
#define SEED time(NULL)

int main(int argc, char **argv)
{
  size_t b = 17;
  size_t n = 22;
  size_t C = 4;
  double alpha[] = { 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 
    0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89};
  size_t N = (1 << n);
  size_t B = (1 << b);
  double error;
  double biterror;
  double n_loops;
  double success;

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

  for (k = 0 ; k < ALEN ; k++)
  {
    /* init trackers */
    n_loops = 0;
    success = 0;
    error = 0;
    biterror = 0;

    size_t K = (size_t)pow(2, alpha[k]*n);

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
      int ncoeff = sfht_deterministic(x, (size_t)N, y, s, K, B, C, L, &loops, &unsat);
    
      
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
    printf("alpha=%f b=%d n=%d error=%f biterror=%f loops=%f success=%f)\n", alpha[k], (int)b, (int)n, error/LOOP, biterror/LOOP, n_loops/LOOP, success/LOOP);
    fprintf(outfile, "%f %d %d %f %d %d %d %d)\n", alpha[k], (int)b, (int)n, error, (int)biterror, (int)n_loops, (int)success, (int)LOOP);
        
    free(y);
    free(s);

  }

  fclose(outfile);

  /* free mem */
  free(x);
  free(index);

  return 0;
}
