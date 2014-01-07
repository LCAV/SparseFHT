/*=========================================================
 * FastHadamard.c - Compute the Hadamard transform
 * of input vector X. If X is a matrix, computes the Hadamard
 * transform of every column in X.
 *
 * [Y] = FastHadamard(X) - Compute the Hadamard transform
 * of input vector X. If X is a matrix, computes the Hadamard
 * transform of every column in X.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2009-2010 The MathWorks, Inc.
 *=======================================================*/
/* $Revision: 1.1.6.2 $ $Date: 2011/01/28 18:11:58 $ */

#include "mex.h"

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "../fht.h"
#include "../sfht.h"

int is_power_of_2(int n)
{
  return !(n & (n-1));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *N, *K, *B, *C, *L, *R, *SEED;  /* pointers to input matrices */
    double *Tfht, *Tsfht;      /* pointers to output matrices */
    int i, j, p;
    size_t *k, *b;
    size_t n, c, l;
    size_t k_len, b_len, r_len, c_len;
    size_t loop, warm, body;
    int max_mag;

    /* check arguments */
    if (nlhs != 2 || nrhs != 7)
      mexErrMsgTxt("HadamardBenchmark: syntax: [Tfht Tsfht] = HadamardBenchmark(N, K, B, C, L, R, SEED)");

    /* get input data */
    N = mxGetPr(prhs[0]); /* first input matrix */
    K = mxGetPr(prhs[1]); /* K: the sparsity */
    k_len = mxGetM(prhs[1])*mxGetN(prhs[1]);
    B = mxGetPr(prhs[2]); /* ALPHA: the over-binning factor */
    b_len = mxGetM(prhs[2])*mxGetN(prhs[2]);
    C = mxGetPr(prhs[3]); /* C: the oversampling factor */
    c_len = mxGetM(prhs[3])*mxGetN(prhs[3]);
    L = mxGetPr(prhs[4]); /* L: the number of iterations of the decoder */
    R = mxGetPr(prhs[5]); /* R: array with sim repetition parameters */
    r_len = mxGetM(prhs[5])*mxGetN(prhs[5]);
    SEED = mxGetPr(prhs[6]); /* SEED: the random number generator seed */

    /* run checks on input */
    if (k_len != b_len)
      mexErrMsgTxt("HadamardBenchmark: syntax: K and B array needs to be the same length.");
    if (k_len != c_len && c_len != 1)
      mexErrMsgTxt("HadamardBenchmark: syntax: C must either be a scalar or the same size as K.");
    if (r_len != 4)
      mexErrMsgTxt("HadamardBenchmark: syntax: R needs to be a length-4 array.");

    /* init RNG */
    srandom((unsigned)SEED[0]);

    /* transform parameters */
    n = (size_t)N[0];
    l = (size_t)L[0];
    k = (size_t *)malloc(k_len*sizeof(size_t));
    b = (size_t *)malloc(b_len*sizeof(size_t));
    for (p = 0 ; p < k_len ; p++)
    {
      k[p] = (size_t)K[p];
      b[p] = (size_t)B[p];
    }

    /* simulation parameters */
    loop = (size_t)R[0];
    warm = (size_t)R[1];
    body = (size_t)R[2];
    max_mag = (int)R[3];

    /* check dimension is a power of two */
    if (!is_power_of_2(n))
      mexErrMsgTxt("SparseFHT2: Transform size needs to be a power of two.");

    /* create output matrix Y */
    plhs[0] = mxCreateDoubleMatrix(loop, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(loop, k_len, mxREAL);
    Tfht = mxGetPr(plhs[0]);
    Tsfht = mxGetPr(plhs[1]);

    /* allocate permutation vector */
    double *x_fht = (double *)malloc(n*sizeof(double));
    for (p = 0 ; p < n ; p++)
      x_fht[p] = (random() % max_mag) * (1 - 2*random()%2);
    double *x = (double *)malloc(n*sizeof(double));
    int *I = (int *)malloc(n*sizeof(int));
    for (i = 0 ; i < n ; i++)
      I[i] = i;

    /*******/
    /* FHT */
    /*******/
    for (j = 0 ; j < loop ; j++)  /* loop for repetition */
    {
      /* now is the time */
      struct timeval t1, t2;
      double elapsedTime;
        
      for (p = 0 ; p < warm ; p++)
        fht(x_fht, n);

      gettimeofday(&t1, NULL);  // start of time interval

      for (p = 0 ; p < body ; p++)
        fht(x_fht, n);

      gettimeofday(&t2, NULL);  // end of time interval

      // compute elapsed time and save it
      elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
      elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
      Tfht[j] = elapsedTime/(double)body;
    }


    /********/
    /* sFHT */
    /********/
    for (i = 0 ; i < k_len ; i++) /* loop for sparsity */
    {
      /* create the 'dummy' output arrays */
      double *y = (double *)malloc(k[i]*sizeof(double));
      int *s = (int *)malloc(k[i]*sizeof(double));

      /* pick right value of C */
      if (c_len == 1)
        c = (size_t)C[0];
      else
        c = (size_t)C[i];

      /* loop for repetition */
      for (j = 0 ; j < loop ; j++)
      {
        /* create new signal vector */
        for (p = 0 ; p < k[i] ; p++)
        {
          /* shuffle support */
          int rn = (random() % n-p) + p;
          int t = I[rn];
          I[rn] = I[p];
          I[p] = t;
          x[I[p]] = (random() % max_mag) * (1 - 2*random()%2);
        }

        /* create its time-domain version */
        fht(x, n);

        /* now is the time */
        struct timeval t1, t2;
        double elapsedTime;
        
        for (p = 0 ; p < warm ; p++)
          sfht_det_opt(x, n, y, s, k[i], b[i], c, l, NULL, NULL);

        gettimeofday(&t1, NULL);  // start of time interval

        for (p = 0 ; p < body ; p++)
          sfht_det_opt(x, n, y, s, k[i], b[i], c, l, NULL, NULL);

        gettimeofday(&t2, NULL);  // end of time interval

        // comput elapsed time and save it
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
        Tsfht[i*loop + j] = elapsedTime / (double)body;;
        
        /* clear signal vector */
        for (p = 0 ; p < k[i] ; p++)
          x[I[p]] = 0;
      }

      /* free memory */
      free(y);
      free(s);
    }

    /* free more memory */
    free(x);
    free(x_fht);
    free(I);

}

