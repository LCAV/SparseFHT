

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "Python.h"
#include "numpy/arrayobject.h"

#include "../../C/fht.h"
#include "../../C/sfht.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#define ALGO_RANDOM 1
#define ALGO_DETERMINISTIC 2
#define ALGO_OPTIMIZED 3

int is_power_of_2(int n)
{
  return !(n & (n-1));
}

/*
 * A python wrapper to call MCMC
 */

  static PyObject*
fht_wrapper (PyObject *dummy, PyObject *args)
{
  PyObject *in=NULL, *out=NULL, *ret=NULL;
  PyArrayObject *iarr=NULL, *oarr=NULL;
  int iarr_nd, oarr_nd;
  npy_intp *iarr_shape, *oarr_shape;
  double *X = NULL, *Y = NULL; /* pointers to input and output matrices */
  int rows, cols;
  int i;

  // Parse the input arguments of the function
  if (!PyArg_ParseTuple(args, "OO!", 
        &in, 
        &PyArray_Type, &out)) return NULL;

  // First argument is the input array
  iarr = (PyArrayObject*)PyArray_FROM_OTF(in, NPY_DOUBLE, NPY_IN_ARRAY);
  if (iarr == NULL) 
  {
    fprintf(stderr, "Could not get pointer to input ndarray");
    goto fail;
  }

  // Second argument is the output array
  oarr = (PyArrayObject*)PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_INOUT_ARRAY);
  if (oarr == NULL)
  {
    fprintf(stderr, "Could not get pointer to output array.");
    goto fail;
  }

  /*vv* code that makes use of arguments *vv*/

  // We run some check on input array
  iarr_nd = PyArray_NDIM(iarr);       // check 2 dimensional
  if (iarr_nd > 2 || iarr_nd == 0)
  {
    fprintf(stderr, "Input array should be 1D or 2D");
    goto fail;
  }
  iarr_shape = PyArray_DIMS(iarr);    // npy_intp array of length nd showing length in each dim.
  if (iarr_shape[0] == 0 || (iarr_nd == 2 && iarr_shape[1] == 0))
  {
    fprintf(stderr, "One of the input dimensions is zero.");
    goto fail;
  }

  // Now some checks on the output array
  oarr_nd = PyArray_NDIM(oarr);
  if (oarr_nd != iarr_nd) 
  {
    fprintf(stderr, "Input and output should have same dimension.");
    goto fail;
  }
  oarr_shape = PyArray_DIMS(oarr);
  if (oarr_shape[0] != iarr_shape[0] || (iarr_nd == 2 && oarr_shape[1] != iarr_shape[1]))
  {
    fprintf(stderr, "Input and output should have same dimensions.");
    goto fail;
  }

  // The transform will be run on the second dimension of the array
  if (iarr_nd == 2)
  {
    rows = iarr_shape[0];
    cols = iarr_shape[1];
  }
  else
  {
    rows = 1;
    cols = iarr_shape[0];
  }

  // recover the data and graph size (finally!)
  X = (double *)PyArray_DATA(iarr);
  Y = (double *)PyArray_DATA(oarr);

  /*^^* code that makes use of arguments *^^*/

  // Copy the data from input to output since the transform is in-place
  for (i = 0 ; i < rows*cols ; i++)
    Y[i] = X[i];

  // Run the FHT
  for (i = 0 ; i < rows ; i++)
    fht(Y + i*cols, cols);

  Py_DECREF(iarr);
  Py_DECREF(oarr);

  // return the energy
  ret = Py_BuildValue("O", out);
  return ret;

fail:
  fprintf(stderr, "Exiting with error.");

  Py_XDECREF(iarr);
  PyArray_XDECREF_ERR(oarr);
  return NULL;
}

  static PyObject*
sfht_wrapper (PyObject *dummy, PyObject *args)
{
  PyObject *arr1=NULL, *arr2=NULL, *arr3=NULL, *arr4=NULL, *arr5=NULL, *ret=NULL;
  PyArrayObject *py_X=NULL, *py_Y=NULL, *py_sup=NULL, *py_U=NULL, *py_I=NULL;
  int X_nd, Y_nd, sup_nd;
  npy_intp *X_shape, *Y_shape, *sup_shape, *U_shape, *I_shape;

  double *X = NULL, *Y = NULL; /* pointers to input and output matrices */
  int *support;                        /* pointer to support array */
  int *U = NULL, *I = NULL; /* pointers to array for performance measure */

  /* parameters */
  int K, B, C, L;
  int T;
  unsigned long Z;
  int n_transform, length;

  /* temporary storage for loops and unsat */
  int loop_num, unsat_chk;
  int *ln=NULL, *us=NULL;

  int i;

  // Parse the input arguments of the function
  if (!PyArg_ParseTuple(args, "OO!O!O!O!iiiiik", 
        &arr1, 
        &PyArray_Type, &arr2,
        &PyArray_Type, &arr3,
        &PyArray_Type, &arr4,
        &PyArray_Type, &arr5,
        &K, &B, &C, &L, &T, &Z)) return NULL;

  /* init RNG for random algorithm */
  if (T == ALGO_RANDOM)
    srandom((unsigned)Z);

  // First argument is the input array
  py_X = (PyArrayObject*)PyArray_FROM_OTF(arr1, NPY_DOUBLE, NPY_IN_ARRAY);
  if (py_X == NULL) 
  {
    fprintf(stderr, "Could not get pointer to input ndarray");
    goto fail;
  }

  // Second argument is the output array
  py_Y = (PyArrayObject*)PyArray_FROM_OTF(arr2, NPY_DOUBLE, NPY_INOUT_ARRAY);
  if (py_Y == NULL)
  {
    fprintf(stderr, "Could not get pointer to output array.");
    goto fail;
  }

  // Third argument is the support array
  py_sup = (PyArrayObject*)PyArray_FROM_OTF(arr3, NPY_INT, NPY_INOUT_ARRAY);
  if (py_sup == NULL)
  {
    fprintf(stderr, "Could not get pointer to output array.");
    goto fail;
  }

  // Fourth argument is the unsat nodes array
  py_U = (PyArrayObject*)PyArray_FROM_OTF(arr4, NPY_INT, NPY_INOUT_ARRAY);
  if (py_U == NULL)
  {
    fprintf(stderr, "Could not get pointer to output array.");
    goto fail;
  }

  // Third argument is the loop number array
  py_I = (PyArrayObject*)PyArray_FROM_OTF(arr5, NPY_INT, NPY_INOUT_ARRAY);
  if (py_I == NULL)
  {
    fprintf(stderr, "Could not get pointer to output array.");
    goto fail;
  }

  // We run some check on input array
  X_nd = PyArray_NDIM(py_X);   // check 2 dimensional
  if (X_nd > 2 || X_nd == 0)
  {
    fprintf(stderr, "Input array should be 1D or 2D");
    goto fail;
  }
  X_shape = PyArray_DIMS(py_X);  // npy_intp array of length nd showing length in each dim.
  if (X_shape[0] == 0 || (X_nd == 2 && X_shape[1] == 0))
  {
    fprintf(stderr, "One of the input dimensions is zero.");
    goto fail;
  }

  // Now some checks on the output array
  Y_nd = PyArray_NDIM(py_Y);
  if (Y_nd != X_nd)
  {
    fprintf(stderr, "Input and output should have same number of dimensions.");
    goto fail;
  }
  Y_shape = PyArray_DIMS(py_Y);
  if ((Y_nd == 2 && (Y_shape[0] != X_shape[0] || Y_shape[1] != K)) || (Y_nd == 1 && Y_shape[0] != K))
  {
    fprintf(stderr, "Mismatch in output dimensions.");
    goto fail;
  }

  // Now some checks on the output array
  sup_nd = PyArray_NDIM(py_sup);
  if (sup_nd != X_nd)
  {
    fprintf(stderr, "Input and output should have same number of dimensions.");
    goto fail;
  }
  sup_shape = PyArray_DIMS(py_sup);
  if (sup_shape[0] != Y_shape[0] || (Y_nd == 2 && sup_shape[1] != Y_shape[1]))
  {
    fprintf(stderr, "Mismatch in support dimensions.");
    goto fail;
  }

  U_shape = PyArray_DIMS(py_U);
  I_shape = PyArray_DIMS(py_I);

  // Some more arguments need to be infered
  if (X_nd == 1)
  {
    n_transform = 1;
    length = X_shape[0];
  }
  else
  {
    n_transform = X_shape[0];
    length = X_shape[1];
  }

  // recover the arrays (finally!)
  X = (double *)PyArray_DATA(py_X);
  Y = (double *)PyArray_DATA(py_Y);
  support = (int32_t *)PyArray_DATA(py_sup);

  if (U_shape[0] == n_transform)
  {
    U = (int32_t *)PyArray_DATA(py_U);
    us = &unsat_chk;
  }
  if (I_shape[0] == n_transform)
  {
    I = (int32_t *)PyArray_DATA(py_I);
    ln = &loop_num;
  }

  /*^^* code that makes use of arguments *^^*/

  /* compute the transform */
  for (i = 0 ; i < n_transform ; i++)  /* over columns if matrix */
  {
    /* compute sparse hadamard transform */
    switch (T)
    {
      case ALGO_RANDOM:
        sfht_random(X+i*length, length, Y+i*K, support+i*K, K, B, C, L,  ln, us);
        break;

      case ALGO_DETERMINISTIC:
        sfht_deterministic(X+i*length, length, Y+i*K, support+i*K, K, B, C, L, ln, us);
        break;

      case ALGO_OPTIMIZED:
      default:
        sfht_det_opt(X+i*length, length, Y+i*K, support+i*K, K, B, C, L, ln, us);
        break;
    }

    /* save performance if required */
    if (us != NULL)
      U[i] = unsat_chk; // number of unsatisfied check nodes
    if (ln != NULL)
      I[i] = loop_num;    // number of loops of decoder

  }

  Py_DECREF(py_X);
  Py_DECREF(py_Y);
  Py_DECREF(py_sup);
  Py_DECREF(py_U);
  Py_DECREF(py_I);

  // return the energy
  ret = Py_BuildValue("");
  return ret;

fail:
  fprintf(stderr, "Exiting with error.");

  Py_XDECREF(py_X);
  PyArray_XDECREF_ERR(py_Y);
  PyArray_XDECREF_ERR(py_sup);
  PyArray_XDECREF_ERR(py_U);
  PyArray_XDECREF_ERR(py_I);
  return NULL;
}

  static PyObject*
benchmark_wrapper (PyObject *dummy, PyObject *args)
{
  PyObject *iarr1=NULL, *iarr2=NULL, *iarr3=NULL, *iarr4=NULL, *oarr1=NULL, *oarr2=NULL, *ret=NULL;
  PyArrayObject *py_K=NULL, *py_B=NULL, *py_C=NULL, *py_R=NULL, *py_Tsfht=NULL, *py_Tfht=NULL;
  int K_nd, B_nd, C_nd, R_nd, Tfht_nd, Tsfht_nd;
  npy_intp *K_shape, *B_shape, *C_shape, *R_shape, *Tfht_shape, *Tsfht_shape;
  int *k = NULL, *b = NULL, *C = NULL, *R = NULL; 
  double *Tfht = NULL, *Tsfht = NULL;      /* pointers to output matrices */
  int n, c, l;
  unsigned int rng_seed;

  size_t k_len, b_len, r_len, c_len;
  size_t loop, warm, body;
  int max_mag;

  int i, j, p;

  // Parse the input arguments of the function
  if (!PyArg_ParseTuple(args, "OOOOO!O!iiI", 
        &iarr1,
        &iarr2,
        &iarr3,
        &iarr4,
        &PyArray_Type, &oarr1,
        &PyArray_Type, &oarr2,
        &n, &l, &rng_seed)) return NULL;

  /* init RNG */
  srandom(rng_seed);

  // Get all nd-arrays
  py_K = (PyArrayObject*)PyArray_FROM_OTF(iarr1, NPY_INT32, NPY_INOUT_ARRAY);
  py_B = (PyArrayObject*)PyArray_FROM_OTF(iarr2, NPY_INT32, NPY_INOUT_ARRAY);
  py_C = (PyArrayObject*)PyArray_FROM_OTF(iarr3, NPY_INT32, NPY_INOUT_ARRAY);
  py_R = (PyArrayObject*)PyArray_FROM_OTF(iarr4, NPY_INT32, NPY_INOUT_ARRAY);
  py_Tsfht = (PyArrayObject*)PyArray_FROM_OTF(oarr1, NPY_DOUBLE, NPY_INOUT_ARRAY);
  py_Tfht = (PyArrayObject*)PyArray_FROM_OTF(oarr2, NPY_DOUBLE, NPY_INOUT_ARRAY);
  if (py_K == NULL || py_B == NULL || py_C == NULL || py_R == NULL || oarr1 == NULL || oarr2 == NULL)
  {
    fprintf(stderr, "Could not get pointer to some array.");
    goto fail;
  }

  K_nd = PyArray_NDIM(py_K);
  B_nd = PyArray_NDIM(py_B);
  C_nd = PyArray_NDIM(py_C);
  R_nd = PyArray_NDIM(py_R);
  Tsfht_nd = PyArray_NDIM(py_Tsfht);
  Tfht_nd = PyArray_NDIM(py_Tfht);

  K_shape = PyArray_DIMS(py_K);
  B_shape = PyArray_DIMS(py_B);
  C_shape = PyArray_DIMS(py_C);
  R_shape = PyArray_DIMS(py_R);
  Tsfht_shape = PyArray_DIMS(py_Tsfht);
  Tfht_shape = PyArray_DIMS(py_Tfht);

  k_len = K_shape[0];
  b_len = B_shape[0];
  c_len = C_shape[0];
  r_len = R_shape[0];

  /* run checks on input */
  if (k_len != b_len)
    fprintf(stderr, "HadamardBenchmark: syntax: K and B array needs to be the same length.");
  if (k_len != c_len && c_len != 1)
    fprintf(stderr, "HadamardBenchmark: syntax: C must either be a scalar or the same size as K.");
  if (r_len != 4)
    fprintf(stderr, "HadamardBenchmark: syntax: R needs to be a length-4 array.");

  /* get the data pointers */
  k = (int *)PyArray_DATA(py_K);
  b = (int *)PyArray_DATA(py_B);
  C = (int *)PyArray_DATA(py_C);
  R = (int *)PyArray_DATA(py_R);
  Tsfht = (double *)PyArray_DATA(py_Tsfht);
  Tfht = (double *)PyArray_DATA(py_Tfht);

  /* simulation parameters */
  loop = (size_t)R[0];
  warm = (size_t)R[1];
  body = (size_t)R[2];
  max_mag = (int)R[3];

  /* check output arrays are large enough */
  if (Tsfht_nd != 2 || Tsfht_shape[0] != loop || Tsfht_shape[1] != k_len || Tfht_shape[0] != loop)
  {
    fprintf(stderr, "The Tsfht output array should be loop x k_len size, and Tfht loop long.");
    goto fail;
  }

  /*^^* code that makes use of arguments *^^*/

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
      Tsfht[i + j*k_len] = elapsedTime / (double)body;;

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


  // Clean up and prepare to exit
  Py_DECREF(py_K);
  Py_DECREF(py_B);
  Py_DECREF(py_C);
  Py_DECREF(py_R);
  Py_DECREF(py_Tsfht);
  Py_DECREF(py_Tfht);

  // return the energy
  ret = Py_BuildValue("OO", oarr1, oarr2);
  return ret;

fail:
  fprintf(stderr, "Exiting with error.");

  Py_XDECREF(py_K);
  Py_XDECREF(py_B);
  Py_XDECREF(py_C);
  Py_XDECREF(py_R);
  PyArray_XDECREF_ERR(py_Tsfht);
  PyArray_XDECREF_ERR(py_Tfht);
  return NULL;
}
static struct PyMethodDef methods[] = {
  {"fht", fht_wrapper, METH_VARARGS, "fht: syntax: fht(x,y)"},
  {"sparse_fht", sfht_wrapper, METH_VARARGS, "sparse_fht: sparse_fht(X, Y, support, U, L, K, B, C, L, T, Z)"},
  {"benchmark", benchmark_wrapper, METH_VARARGS, "HadamardBenchmark: benchmark(K, B, C, R, Tsfht, Tfht, N, L, SEED)"},
  {NULL, NULL, 0, NULL}
};

  PyMODINIT_FUNC
initsparsefht_wrapper (void)
{
  PyObject *m;

  m = Py_InitModule("sparsefht_wrapper", methods);
  import_array();

  PyModule_AddObject(m, "ALGO_RANDOM", PyInt_FromLong(ALGO_RANDOM));
  PyModule_AddObject(m, "ALGO_DETERMINISTIC", PyInt_FromLong(ALGO_DETERMINISTIC));
  PyModule_AddObject(m, "ALGO_OPTIMIZED", PyInt_FromLong(ALGO_OPTIMIZED));
}
