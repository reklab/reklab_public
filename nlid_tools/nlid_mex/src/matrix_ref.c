/*

  matrix_ref.c
  This function returns the offset from the base-pointer for the 
  [i,j,k] element in and n by n by n, third-order tensor.

  offset = matrix_ref ( n,i,j,k);

  David Westwick   August 24, 1993
 
*/

#ifdef __STDC__
int matrix_ref (
  int n,					    /* size of matrix */
  int i,				/* indeces of matrix  element */
  int j,
  int k
  )
#else
int matrix_ref (n,i,j,k)
  int n, i, j, k;
#endif
{
  return i * (n * n) +  j * n + k;
  }



