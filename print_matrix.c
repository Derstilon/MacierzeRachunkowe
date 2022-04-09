#include <stdio.h>

#define A(i, j) a[(j)*lda + (i)]

void print_matrix(int m, int n, double **a, int lda)
{
   int i, j;
   for (i = 0; i < m; i++)
   {
      for (j = 0; j < n; j++)
         printf("%le ", a[i][j]);
      printf("\n");
   }
   printf("\n");
}
