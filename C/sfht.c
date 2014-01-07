
#include "sfht.h"
#include "gf2.h"
#include "fht.h"

#define EPSILON_0 1e-5
#define EPSILON_2 1e-8

#include <math.h>

#include <stdlib.h>
#include <time.h>

/*
 * Sparse Fast Hadamard Transform
 * using Random hash matrices
 *
 * Arguments:
 * x: input vector  (size n)
 * n: the transform size (size of x too)
 * y: output vector (size k)
 * s: support vector (size k)
 * k: the sparsity (and size of y)
 * B: the number of bins
 * C: oversampling factor
 * L: number of iterations of decoder
 * loops_done: a pointer to an integer to store the number of loops used by the peeling decoder.
 * unsat: a pointer to an integer to store the number of unsatisfied check nodes.
 *
 * Return: the number of coefficients found.
 */
int sfht_random(double *x, size_t N, double *y, int *s, size_t K, size_t B, int C, int L, int *loops_done, int *unsat)
{
  int n_coeff_found = 0;
  int c;

  /* the randomization parameters */
  gf2_matrix_t *U[C], *V[C], *Vinv[C];

  /* set all parameters */
  size_t n = (size_t)log2(N);

  /* allocate arrays */
  for (c = 0 ; c < C ; c++)
  {
    U[c] = gf2_alloc_mat(n);
    V[c] = gf2_alloc_mat(n);
    Vinv[c] = gf2_alloc_mat(n);

    /* matrix randomness */
    gf2_rand_mat(U[c], Vinv[c]);
    gf2_copy(U[c], V[c]);
    gf2_transpose(U[c]);
  }

  /* run SFHT algorithm */
  n_coeff_found = sfht_base(x, N, y, s, K, B, C, U, V, Vinv, L, loops_done, unsat);

  /* free allocated stuff */
  for (c = 0 ; c < C ; c++)
  {
    gf2_free_mat(U[c]);
    gf2_free_mat(V[c]);
    gf2_free_mat(Vinv[c]);
  }

  /* return */
  return n_coeff_found;

}

/*
 * Sparse Fast Hadamard Transform
 * using deterministic sampling pattern
 *
 * Arguments:
 * x: input vector  (size n)
 * N: the transform size (size of x too)
 * y: output vector (size k)
 * s: support vector (size k)
 * k: the sparsity (and size of y)
 * B: number of buckets
 * C: oversampling factor
 * L: maximum number of iterations of decoder
 * l: a parameter for the deterministic hash pattern
 * hash_type: choose the hashing scheme used
 * loops_done: a pointer to an int where the number of loops run is stored
 * unsat_chk: a pointer to an int to store the number of unsatisfied check nodes not used if null
 *
 * Return: the number of unsatisfied checks
 */
int sfht_deterministic(double *x, size_t N, double *y, int *s, size_t K, size_t B, 
    int C, int L, int *loops_done, int *unsat)
{

  int c, i;
  int n_coeff_found = 0;
  size_t a;

  /* the randomization parameters */
  //gf2_matrix_t *U[C], *V[C], *Vinv[C];
  gf2_matrix_t **U, **V, **Vinv;
  U = (gf2_matrix_t **)malloc(C*sizeof(gf2_matrix_t *));
  V = (gf2_matrix_t **)malloc(C*sizeof(gf2_matrix_t *));
  Vinv = (gf2_matrix_t **)malloc(C*sizeof(gf2_matrix_t *));

  /* set all parameters */
  size_t n = (size_t)log2(N);
  size_t b = (size_t)log2(B);

  /* allocate arrays */
  for (c = 0 ; c < C ; c++)
  {
    U[c] = gf2_alloc_mat(n);
    V[c] = gf2_alloc_mat(n);
    Vinv[c] = gf2_alloc_mat(n);

    /* pick deterministic matrix */
    a = (size_t)(c*n/(double)C+0.5);
    
    for (i = 0 ; i < n-b ; i++)
      V[c]->rows[i] = (1 << ((a+b+i)%n));
    for (i = n-b ; i < n ; i++)
      V[c]->rows[i] = (1 << ((i-n+b+a)%n));

    gf2_copy(V[c], U[c]);
    gf2_transpose(U[c]);
    gf2_copy(U[c], Vinv[c]);
  }
  //mexPrintf("\n");

  /* run SFHT algorithm */
  n_coeff_found = sfht_base(x, N, y, s, K, B, C, U, V, Vinv, L, loops_done, unsat);

  /* free allocated stuff */
  for (c = 0 ; c < C ; c++)
  {
    gf2_free_mat(U[c]);
    gf2_free_mat(V[c]);
    gf2_free_mat(Vinv[c]);
  }
  free(U);
  free(V);
  free(Vinv);

  return n_coeff_found;

}

/*
 * Sparse Fast Hadamard Transform
 * Algorithm 2 with 'OFDM trick'
 *
 * Arguments:
 * x: input vector  (size n)
 * n: the transform size (size of x too)
 * y: output vector (size k)
 * s: support vector (size k)
 * k: the sparsity (and size of y)
 * B: the number of bins
 * C: oversampling factor
 * U: The C random matrices in F_2^{n x n}
 * V: Their transposes V=U^T
 * Vinv: The inverses of the transposes Uinv = U^{-T}
 * L: maximum number of iterations of decoder
 * loops_done: a pointer to an integer to store the number of loops used by the peeling decoder.
 * unsat: a pointer to an integer to store the number of unsatisfied check nodes.
 *
 * Return: the number of coefficients found.
 */
int sfht_base(double *x, size_t N, double *y, int *s, size_t K, 
    size_t B, int C, gf2_matrix_t **U, gf2_matrix_t **V, gf2_matrix_t **Vinv, int L, int *loops_done, int *unsat)
{
  int i, j, l, c, d;

  /* set all parameters */
  size_t n = (size_t)log2(N);
  size_t b = (size_t)log2(B);
  size_t D = n-b+1;
  double normalization = sqrt(N)/B;

  /* sample */

  /* vector to store all the hashes */
  double *hash = (double *)malloc(C*D*B*sizeof(double));

  /* compute the hashes, also with modulation */
  for (c = 0 ; c < C ; c++)
  {
    double *h = hash + c*D*B;

    fast_hadamard_hash(x, N, h, B, U[c], 0);
    for (j = 0 ; j < D-1 ; j++)
      fast_hadamard_hash(x, N, h+(j+1)*B, B, U[c], V[c]->rows[j]);

  }

  /* Peeling decoder */

  /* initialize number of recovered coefficients */
  int k = 0;

  /* main loop */
  for (l = 0 ; l < L ; l++)
  {
    /* check when all is done */
    int singleton_found = 0;

    for (c = 0 ; c < C ; c++)
    {
      /* get pointer to start of this part of tree */
      double *h = hash + c*B*D;

      for (i = 0 ; i < B ; i++)
      {
        /* skip if bucket is empty */
        for (j = 0 ; j < D ; j++)
          if (fabs(h[i+j*B]) > EPSILON_0)
            break;
        if (j == D)
          continue;
        //if (fabs(h[i]) < EPSILON_0)
          //continue;

        /* recover the index bit by bit (a.k.a OFDM trick) */
        gf2_vector_t index = 0;
        int success = 1;
        for (j = 0 ; j < D-1 ; j++)
        {
          double r = h[(j+1)*B+i]/h[i];
          if (fabs(fabs(r)-1) < EPSILON_2)
          {
            if (r < 0)
              index ^= (1 << j);
          }
          else
          {
            success = 0;
            break;
          }
        }

        /* skip if more than one coeff in bucket */
        if (success == 0)
          continue;

        /* a singleton was found */
        singleton_found = 1;

        /* store the newly found coefficient */
        y[k] = normalization*h[i];
        s[k] = gf2_mat_vec(Vinv[c], (index) ^ (i << (n-b)));

        /* peel the recovered coefficient from other check nodes */
        for (j = 0 ; j < c ; j++)
        {
          double *h2 = hash + j*D*B;
          gf2_vector_t v = gf2_mat_vec(V[j], s[k]) >> (n-b);
          h2[v] -= h[i];
          for (d = 1 ; d < D ; d++)
            h2[d*B+v] -= h[i]*(1 - 2*gf2_parity(V[j]->rows[d-1] & s[k]));
        }

        for (j = c+1 ; j < C ; j++)
        {
          double *h2 = hash + j*D*B;
          gf2_vector_t v = gf2_mat_vec(V[j], s[k]) >> (n-b);
          h2[v] -= h[i];
          for (d = 1 ; d < D ; d++)
            h2[d*B+v] -= h[i]*(1 - 2*gf2_parity(V[j]->rows[d-1] & s[k]));
        }

        /* set current element to zero */
        h[i] = 0.;

        /* increment k */
        k++;

        /* stopping criterion */
        if (k == K)
          goto end;
      }
    }

    /* exit if there is no singleton node left */
    if (!singleton_found)
      goto end;
  }

end:

  /* count unsat check nodes */
  if (unsat != NULL)
  {
    *unsat = 0;
    for (c = 0 ; c < C ; c++)
      for (i = 0 ; i < B ; i++)
        if (fabs(hash[c*D*B+i]) >= EPSILON_0)
          (*unsat)++;
  }

  /* number of loops */
  if (loops_done != NULL)
    *loops_done = l+1;

  free(hash);

  /* return number of coefficients found */
  return k;

}

/* Sparse Hadamard hashing function */
/* the output needs normalized by sqrt(n)/B */
void fast_hadamard_hash(double *x, size_t n, double *hash, size_t B, gf2_matrix_t *sigma, gf2_vector_t a)
{
  int i;
  int bin_size = n/B;

  /* permute and downsample data in time-domain */
  for (i = 0 ; i < B ; i++)
  {
    gf2_vector_t t = gf2_mat_vec(sigma, i*bin_size);
    hash[i] = x[t ^ a];
  }

  /* Hadamard transform */
  fht(hash, B);

}

/*
 * Sparse Fast Hadamard Transform
 * using deterministic sampling pattern
 * and optimized bitwise operations to perform them
 *
 * Arguments:
 * x: input vector  (size n)
 * N: the transform size (size of x too)
 * y: output vector (size k)
 * s: support vector (size k)
 * k: the sparsity (and size of y)
 * B: number of buckets
 * C: oversampling factor
 * L: maximum number of iterations of decoder
 * l: a parameter for the deterministic hash pattern
 * hash_type: choose the hashing scheme used
 * loops_done: a pointer to an int where the number of loops run is stored
 * unsat_chk: a pointer to an int to store the number of unsatisfied check nodes not used if null
 *
 * Return: the number of unsatisfied checks
 */
int sfht_det_opt(double *x, size_t N, double *y, int *s, size_t K, size_t B, int C, int L, int *loops_done, int *unsat)
{

  int i, j, c, d;
  uint32_t a;
  int loop;

  /* set all parameters */
  size_t n = (size_t)log2(N);
  size_t b = (size_t)log2(B);
  size_t D = n-b+1;
  double normalization = sqrt(N)/B;

  //size_t cond = n/2;
  //size_t cond = (C-1)*n/C;

  /* sample */

  /* vector to store all the hashes */
  double *hash = (double *)malloc(C*D*B*sizeof(double));
  if (hash == NULL)
  {
    fprintf(stderr, "SFHT_det: Error could not allocate memory.\n");
    exit(1);
  }

  /* compute the hashes, also with modulation */
  for (c = 0 ; c < C ; c++)
  {
    double *h = hash + c*D*B;

    //a = (b <= cond) ? c*n/C : (c+1)*(n-b);
    a = select_hash(b, n, c, C);

    /* Hash without modulation */
    for (i = 0 ; i < B ; i++)
      h[i] = x[inv_hashing(i, a, b, n)]; /* downsample and permute */
    fht(h, B);                       /* Hadamard transform */

    /* Hash with modulation */
    for (j = 0 ; j < D-1 ; j++)
    {
      for (i = 0 ; i < B ; i++)
        h[(j+1)*B+i] = x[inv_hashing(i, a, b, n) ^ (1 << null_vector(j, a, b, n))]; /* downsample, permute, delay */
      fht(h+(j+1)*B, B);                                             /* Hadamard transform */
    }
  }

  /* Peeling decoder */

  /* initialize number of recovered coefficients */
  int k = 0;

  /* main loop */
  for (loop = 0 ; loop < L ; loop++)
  {
    /* check when all is done */
    int singleton_found = 0;

    /* loop over the C hash */
    for (c = 0 ; c < C ; c++)
    {
      /* get pointer to start of this part of graph */
      double *h = hash + c*B*D;

      /* the hash factor */
      //a = (b <= cond) ? c*n/C : (c+1)*(n-b);
      a = select_hash(b, n, c, C);

      for (i = 0 ; i < B ; i++)
      {
        /* skip if bucket is empty */
        if (fabs(h[i]) < EPSILON_0)
          continue;

        /* recover the index bit by bit (a.k.a OFDM trick) */
        uint32_t index = inv_hashing(i, a, b, n);
        int success = 1;
        for (j = 0 ; j < D-1 ; j++)
        {
          double r = h[(j+1)*B+i]/h[i];
          if (fabs(fabs(r)-1) < EPSILON_2)
          {
            if (r < 0)
              index ^= (1 << null_vector(j, a, b, n));
          }
          else
          {
            success = 0;
            break;
          }
        }

        /* skip if more than one coeff in bucket */
        if (success == 0)
          continue;

        /* a singleton was found */
        singleton_found = 1;

        /* store the newly found coefficient */
        y[k] = normalization*h[i];
        s[k] = index;

        /* peel the recovered coefficient from other check nodes */
        for (j = 0 ; j < c ; j++)
        {
          double *h2 = hash + j*D*B;
          //uint32_t a2 = (b <= cond) ? j*n/C : (j+1)*(n-b);
          uint32_t a2 = select_hash(b, n, j, C);
          uint32_t v = hashing(s[k], a2, b, n);
          h2[v] -= h[i];
          for (d = 1 ; d < D ; d++)
            h2[d*B+v] -= h[i]*(1 - 2*((s[k] >> null_vector(d-1, a2, b, n)) & 1));
        }

        for (j = c+1 ; j < C ; j++)
        {
          double *h2 = hash + j*D*B;
          //uint32_t a2 = (b <= cond) ? j*n/C : (j+1)*(n-b);
          uint32_t a2 = select_hash(b, n, j, C);
          uint32_t v = hashing(s[k], a2, b, n);
          h2[v] -= h[i];
          for (d = 1 ; d < D ; d++)
            h2[d*B+v] -= h[i]*(1 - 2*((s[k] >> null_vector(d-1, a2, b, n)) & 1));
        }

        /* set current element to zero */
        h[i] = 0.;

        /* increment k */
        k++;

        /* stopping criterion */
        if (k == K)
          goto end;
      }
    }

    /* exit if there is no singleton node left */
    if (!singleton_found)
      goto end;
  }

end:

  /* count unsat check nodes */
  if (unsat != NULL)
  {
    *unsat = 0;
    for (c = 0 ; c < C ; c++)
      for (i = 0 ; i < B ; i++)
        if (fabs(hash[c*D*B+i]) >= EPSILON_0)
          (*unsat)++;
  }

  /* free the hash */
  free(hash);

  /* number of loops */
  if (loops_done != NULL)
    *loops_done = loop+1;

  return k;

}

/*
 * Low regime
 * 0 < \alpha < (C-1)/C
 */
uint32_t select_hash(uint32_t b, uint32_t n, uint32_t c, uint32_t C)
{
  return (uint32_t)(c*n/(double)C+0.5);
}

uint32_t hashing(uint32_t val, uint32_t a, uint32_t b, uint32_t n)
{
  return ((val >> a) | (val << (n-a))) & ((1 << b)-1);
}

uint32_t inv_hashing(uint32_t val, uint32_t a, uint32_t b, uint32_t n)
{
  return ((val << a) | (val >> (n-a))) & ((1 << n)-1);
}

uint32_t null_vector(uint32_t d, uint32_t a, uint32_t b, uint32_t n)
{
  return (a+b+d) % n;
}

