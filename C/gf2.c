
#include "gf2.h"

/* look-up table for parity check */
static const uint8_t ParityTable256[256] = 
{
#define P2(n) n, n^1, n^1, n
#define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
  P6(0), P6(1), P6(1), P6(0)
};

#if 0
gf2_vector_t (*gf2_parity_func[4])(gf2_vector_t) = { &gf2_parity_8, &gf2_parity_16, &gf2_parity_32, &gf2_parity_32 };
#endif

/* 
 * Create a random non-singular matrix U of size NxN in GF(2) and its inverse V.
 */ 
gf2_matrix_t *gf2_rand_mat(gf2_matrix_t *U, gf2_matrix_t *V)
{
  int i,j;

  /* check size */
  if (U->size != V->size)
  {
    fprintf(stderr, "gf2_rand_mat: Error, sizes are different.\n");
    return NULL;
  }

  /* fill the matrix with random bits */
  while (1)
  {
    for (i = 0 ; i < U->size ; i++)
    {
      /* draw a random number between 1 and 2^N-1 without repetition */
      /* XXX The algo is ineficient, but N is always small XXX */
      while (1)
      {
        U->rows[i] = (gf2_vector_t)(random() % ((1ul << U->size)-1)) + 1;
        /* check all previously drawn elements */
        for (j = 0 ; j < i ; j++)
          if (U->rows[i] == U->rows[j])
            break;
        /* if reach current element without finding repetition, exit loop */
        if (j == i)
          break;
      }
    }
    /* repeat the loop, unless the matrix is full rank */
    if (gf2_inverse(U, V) != NULL)
      break;
  }

  return U;
}

void gf2_rand_nonsing(gf2_matrix_t *A)
{
  int a, t;

  /* create new matrix T */
  gf2_matrix_t *T = gf2_alloc_mat(A->size);
  T->size = A->size;

  /* assign matrices to zero */
  for (t = 0 ; t < T->size ; t++)
  {
    T->rows[t] = 0;
    A->rows[t] = 0;
  }

  /* follow Dana Randall's algorithm to assign A and T in a non-singular fashion */
  for (a = 0 ; a < A->size ; a++)
  {
    /* pick random line */
    gf2_vector_t v = (random() % ((1 << (A->size - a)) - 1)) + 1;
    int r = -1;

    /* look for current line and also set bits in T's row */
    for (t = 0 ; t < T->size ; t++)
    {
      /* skip rows already occupied */
      if (T->rows[t] != 0)
        continue;

      /* we just reached the first one in v */
      if (v & 1)
      {
        if (r < 0)
        {
          r = t;
          //A->rows[a] = (1 << t); /* echelon row to add in A */
          A->rows[a] = t;
        }
        /* finish assigning the non-zero bits */
        T->rows[r] |= 1 << t;
      }
      v >>= 1;
    }
  }

  for (a = 0 ; a < A->size ; a++)
    A->rows[a] = T->rows[A->rows[a]];

  /* free matrix T */
  gf2_free_mat(T);
}

/* allocate GF(2) matrix of size NxN */
gf2_matrix_t *gf2_alloc_mat(size_t N)
{
  /* only allows matrix up to 32 by 32 */
  if (N > GF2_MAX_SIZE)
  {
    fprintf(stderr, "gf2_alloc_mat: Error, matrix larger than %dx%d.\n", GF2_MAX_SIZE, GF2_MAX_SIZE);
    return NULL;
  }

  /* allocate structure */
  gf2_matrix_t *U = (gf2_matrix_t *)malloc(sizeof(gf2_matrix_t));
  if (U == NULL)
  {
    fprintf(stderr, "gf2_alloc_mat: Error, could not allocate memory.\n");
    return NULL;
  }

  /* allocate rows */
  U->rows = (gf2_vector_t *)malloc(N*sizeof(gf2_vector_t));
  if (U->rows == NULL)
  {
    fprintf(stderr, "gf2_alloc_mat: Error, could not allocate memory.\n");
    free(U);
    return NULL;
  }

  U->size = N;

  return U;
}

/* free the memory allocated to U */
void gf2_free_mat(gf2_matrix_t *U)
{
  if (U != NULL && U->rows != NULL)
    free(U->rows);
  if (U != NULL)
    free(U);
}

/* 
 * computes the inverse of U in GF(2) and stores it in V
 * - assumes U is full-rank
 */
gf2_matrix_t *gf2_inverse(gf2_matrix_t *U, gf2_matrix_t *V)
{
  int r, c;

  /* check size */
  if (U->size != V->size)
  {
    fprintf(stderr, "gf2_inverse: Error, sizes are different.\n");
    return NULL;
  }

  /* allocate matrix */
  gf2_matrix_t *T = gf2_alloc_mat(U->size);

  /* make T a copy of U */
  T->size = U->size;
    for (r = 0 ; r < T->size ; r++)
      T->rows[r] = U->rows[r];

  /* allocate the matrix V as identity of size N */
  V->size = U->size;
  for (r = 0 ; r < V->size ; r++)
    V->rows[r] = 1 << r;

  /* compute the inverse using direct Gauss-Jordan elimination */
  for (c = 0 ; c < T->size ; c++)
  {
    /* find a pivot, if necessary */
    if (((T->rows[c] >> c) & 1) == 0)
    {
      for (r = c ; r < T->size ; r++)
      {
        if (((T->rows[r] >> c) & 1) == 1)
        {
          /* swap row r and c */
          gf2_vector_t tmp = T->rows[c];
          T->rows[c] = T->rows[r];
          T->rows[r] = tmp;
          tmp = V->rows[c];
          V->rows[c] = V->rows[r];
          V->rows[r] = tmp;
          /* break */
          break;
        }
      }
      /* if we reach this, the algorithm failed */
      if (c == T->size)
      {
        gf2_free_mat(T);
        return NULL;
      }
    }

    /* elimination phase */
    for (r = 0 ; r < T->size ; r++)
    {
      if (r != c && ((T->rows[r] >> c) & 1) == 1)
      {
        T->rows[r] ^= T->rows[c];
        V->rows[r] ^= V->rows[c];
      }
    }
  }

  /* check that the matrix is not rank deficient */
  for (r = 0 ; r < T->size ; r++)
  {
    /* if there is a zero line, the matrix is rank deficient */
    if (T->rows[r] == 0)
    {
      gf2_free_mat(T);
      return NULL;
    }
  }

  /* free temporary matrix T */
  gf2_free_mat(T);

  /* return the pointer to inverse matrix */
  return V;
}

/* Transposes U */
void gf2_transpose(gf2_matrix_t *U)
{
  int r, c;

  /* make a copy of U */
  gf2_matrix_t *C = gf2_alloc_mat(U->size);
  C->size = U->size;
  for (r = 0 ; r < U->size ; r++)
  {
    C->rows[r] = U->rows[r];
    U->rows[r] = 0; // zeros U at the same time
  }

  /* transpose U */
  for (r = 0 ; r < U->size ; r++)
    for (c = 0 ; c < U->size ; c++)
    {
      U->rows[c] |= (C->rows[r] & 1) << r;
      C->rows[r] >>= 1;
    }

  /* free the copy */
  gf2_free_mat(C);

}

/* 
 * multiplies an array of size len of vectors V in GF(2) 
 * stored in an integer by matrix U
 * infers the length of V from the size of U 
 */
gf2_vector_t gf2_mat_vec(gf2_matrix_t *U, gf2_vector_t V)
{
  int r;
  gf2_vector_t p = 0;
  //gf2_vector_t (*p_func)(gf2_vector_t) = gf2_parity_func[U->size/8];

  for (r = 0 ; r < U->size ; r++)
  {
    p |= gf2_parity(U->rows[r] & V) << r;
    //p |= ((*p_func)(U->rows[r] & V)) << r;
  }
    
  return p;
}

/*
 * compute the number of ones in the vector V
 */
int gf2_parity(gf2_vector_t V)
{
  /* little more advanced parity check O(number of non-zero bits) */
#if 0
  int parity = 0; /* initialize parity check */
  while (V)
  {
    parity ^= 1;  /* count a bit */
    V &= V-1;     /* this anihilates the left-most non-zero bit */
  }        
  return parity;
#else
  V ^= (V >> 16);
  V ^= (V >> 8);
  return ParityTable256[V & 0xff];
#endif
}

#if 0
gf2_vector_t gf2_parity_32(gf2_vector_t V)
{
  V ^= (V >> 16);
  V ^= (V >> 8);
  return ParityTable256[V & 0xff];
}

gf2_vector_t gf2_parity_16(gf2_vector_t V)
{
  V ^= V >> 8;
  return ParityTable256[V & 0xff];
}

gf2_vector_t gf2_parity_8(gf2_vector_t V)
{
  return ParityTable256[V];
}
#endif

/*
 * copy matrix U into V
 */
void gf2_copy(gf2_matrix_t *U, gf2_matrix_t *V)
{
  int r;
  V->size = U->size;
  for (r = 0 ; r < U->size ; r++)
    V->rows[r] = U->rows[r];
}

/*
 * compute product of matrix and a number of vectors, in-place
 */
void gf2_mat_vec_array(gf2_matrix_t *U, gf2_vector_t *V, size_t len)
{
  while (len-- > 0)
    V[len] = gf2_mat_vec(U, V[len]);
}

/* simple printing routines */
void gf2_print_vec(gf2_vector_t V, size_t N)
{
  printf("[ ");
  while (N > 1)
  {
    printf("%d ", V & 1);
    V >>= 1;
    N--;
  }
  printf("%d ]\n", V & 1);
}

void gf2_print_mat(gf2_matrix_t *U)
{
  int r;
  for (r = 0 ; r < U->size ; r++)
    gf2_print_vec(U->rows[r], U->size);
}

