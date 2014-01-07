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
#include "../fht.h"

int is_power_of_2(int n)
{
  int l = 0;
  int m = n;
  while (m != 0)
  {
    m >>= 1;
    l += 1;
  }
  if (n != (1 << (l-1)))
    return 0;
  else
    return 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *X;     /* pointers to input matrices */
    double *Y;     /* pointers to output matrices */
    int i, rows, cols;

    /* get input data */
    X = mxGetPr(prhs[0]); /* first input matrix */

    /* dimension of input matrix */
    rows = mxGetM(prhs[0]);
    cols = mxGetN(prhs[0]);

    /* check dimension is a power of two */
    if ((rows != 1 && !is_power_of_2(rows)) || (rows == 1 && !is_power_of_2(cols)))
      mexErrMsgTxt("FastHadamard: transform size needs to be a power of two.");

    /* create output matrix Uv */
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    Y = mxGetPr(plhs[0]);
    for (i=0 ; i < rows*cols ; i++)
      Y[i] = X[i];

    /* compute the transform */
    if (rows == 1)  /* over rows if array */
      fht(Y, cols);
    else
      for (i = 0 ; i < cols ; i++)  /* over columns if matrix */
        fht(Y+i*rows, rows);
}
