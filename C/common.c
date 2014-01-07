
#include "common.h"
#include <stdio.h>

/* 
 * Sort a sparse array in increasing order of its indices
 *
 * value is the value array
 * index is the index array
 * l is the index of the first element to sort
 * r is the index of the last element to sort
 */
void quicksort(int *array, int l, int r)
{
  int i,j;
  int key;
  
  if( l < r)
  {

    int tmp = array[r];
    array[r] = array[(l+r)/2];
    array[(l+r)/2] = tmp;

    key = array[r];

    i = l;
    j = r-1;

    while(i <= j)
    {
      while (i <= r && array[i] <= key)
      {
        i++;
        continue;
      }

      while (j >= l && array[j] > key)
      {
        j--;
        continue;
      }

      /* swap two elements */
      if (i < j)
      {
        int t_ind = array[i];
        array[i] = array[j];
        array[j] = t_ind;
      }
    }

    tmp = array[r];
    array[r] = array[j+1];
    array[j+1] = tmp;

    /* recursively sort the lesser list */
    quicksort(array, l, j);
    quicksort(array, j+2, r);
  }
}

void swap2(double* value, gf2_vector_t *index, int k, int l)
{
  double t1 = value[k];
  value[k] = value[l];
  value[l] = t1;
  gf2_vector_t t2 = index[k];
  index[k] = index[l];
  index[l] = t2;
}

/* 
 * Sort a sparse array in increasing order of its indices
 *
 * value is the value array
 * index is the index array
 * l is the index of the first element to sort
 * r is the index of the last element to sort
 */
void sparse_quicksort(double *value, gf2_vector_t *index, int l, int r)
{
  int i,j;
  gf2_vector_t key;
  
  if( l < r)
  {
    /* choose pivot in middle and swap */
    swap2(value, index, (l+r)/2, r);

    /* set key parameters */
    key = index[r];
    i = l;
    j = r-1;

    while(i <= j)
    {
      while (i <= r && index[i] <= key)
        i++;
      while (j >= l && index[j] > key)
        j--;

      /* swap two elements */
      if (i < j)
        swap2(value, index, i, j);
    }

    /* swap pivot back */
    swap2(value, index, j+1, r);

    /* recursively sort the lesser list */
    sparse_quicksort(value, index, l, j);
    sparse_quicksort(value, index, j+2, r);
  }
}

