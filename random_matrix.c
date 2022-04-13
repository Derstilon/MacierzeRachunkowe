#include <stdlib.h>

void random_matrix(int m, int n, double **a, int lda)
{
   srand48(2);
   double drand48();
   int i, j;

   for (j = 0; j < n; j++)
      for (i = 0; i < m; i++)
         a[i][j] = 2.0 * drand48() - 1.0;
}
