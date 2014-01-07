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
#include "../sfht.h"

#define ALGO_RANDOM 1
#define ALGO_DETERMINISTIC 2
#define ALGO_OPTIMIZED 3

int is_power_of_2(int n)
{
  return !(n & (n-1));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *X;              /* pointers to input matrices */
    double *Y, *S, *U, *I;  /* pointers to output matrices */
    int *support;
    int loop_num, unsat_chk;
    int *ln=NULL, *us=NULL;
    int i, j, rows, cols;
    size_t K, B, C, L, G, H;
    int T;
    unsigned long Z;
    size_t n_transform, length;

    /* check arguments */
    if (nlhs < 2 || nlhs > 4 || nrhs != 7)
      mexErrMsgTxt("SparseFHT_mex: syntax: [Y, S, U, I] = SparseFHT_mex(X, K, B, C, L, T, Z) \n\
\n\
Sparse Fast Hadamard Transform\n\
with deterministic sampling pattern\n\
\n\
Input arguments:\n\
X: input vector  (size n)\n\
K: the sparsity (and size of y)\n\
B: number of buckets\n\
C: oversampling factor\n\
L: maximum number of iterations of decoder\n\
T: Type of algorithm to use (random (1) / deterministic (2) / optimized (3))\n\
Z: Seed for the RNG (only used for random algorithm)\n\
\n\
Output arguments: \n\
Y: output vector (size K)\n\
S: support vector (size K)\n\
U: the number of unsatisfied checks (optional)\n\
I: the number of loops run (optional)\n");

    /* get input data */
    X = mxGetPr(prhs[0]);            /* The input signal */
    K = (size_t)mxGetPr(prhs[1])[0]; /* K: the sparsity */
    B = (size_t)mxGetPr(prhs[2])[0]; /* B: the number of buckets */
    C = (size_t)mxGetPr(prhs[3])[0]; /* C: the oversampling factor */
    L = (size_t)mxGetPr(prhs[4])[0]; /* L: the maximum number of iterations of the decoder */
    T = (int)mxGetPr(prhs[5])[0];    /* T: Type of algorithm */
    Z = (unsigned long)mxGetPr(prhs[6])[0]; /* Z: Seed for RNG */

    /* init RNG for random algorithm */
    if (T == ALGO_RANDOM)
      srandom((unsigned)Z);

    /* dimension of input matrix */
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);

    /* check dimension is a power of two */
    if ((rows != 1 && !is_power_of_2(rows)) || (rows == 1 && !is_power_of_2(cols)))
      mexErrMsgTxt("SparseFHTDet: Dimension needs to be a power of two.");

    /* allocate integer support array */
    support = (int *)malloc(K*sizeof(int));
    if (support == NULL)
      mexErrMsgTxt("SparseFHTDet: Not enough memory.");

    /* create output matrix Y */
    if (rows == 1)
    {
      n_transform = 1;
      length = cols;
      plhs[0] = mxCreateDoubleMatrix(n_transform, K, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(n_transform, K, mxREAL);
    }
    else
    {
      n_transform = cols;
      length = rows;
      plhs[0] = mxCreateDoubleMatrix(K, n_transform, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(K, n_transform, mxREAL);
    }
    Y = mxGetPr(plhs[0]);
    S = mxGetPr(plhs[1]);

    /* optional arrays for performance metric */
    if (nlhs >= 3)
    {
      plhs[2] = mxCreateDoubleMatrix(1, n_transform, mxREAL);
      U = mxGetPr(plhs[2]);
      us = &unsat_chk;
    }
    if (nlhs >= 4)
    {
      plhs[3] = mxCreateDoubleMatrix(1, n_transform, mxREAL);
      I = mxGetPr(plhs[3]);
      ln = &loop_num;
    }

    /* compute the transform */
    for (i = 0 ; i < n_transform ; i++)  /* over columns if matrix */
    {
      /* compute sparse hadamard transform */
      switch (T)
      {
        case ALGO_RANDOM:
          sfht_random(X+i*length, length, Y+i*K, support, K, B, C, L, ln, us);
          break;

        case ALGO_DETERMINISTIC:
          sfht_deterministic(X+i*length, length, Y+i*K, support, K, B, C, L, ln, us);
          break;

        case ALGO_OPTIMIZED:
        default:
          sfht_det_opt(X+i*length, length, Y+i*K, support, K, B, C, L, ln, us);
          break;
      }

      /* save performance if required */
      if (us != NULL)
        U[i] = unsat_chk; // number of unsatisfied check nodes
      if (ln != NULL)
        I[i] = loop_num;  // number of loops of decoder

      /* copy support to output array */
      for (j = 0 ; j < K ; j++)
        if (support[j] >= 0 && support[j] < length)
          S[i*K+j] = (double)support[j];
        else
          S[i*K+j] = -1;
    }

    /* free memory */
    free(support);
}
