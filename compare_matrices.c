#define abs(x) ((x) < 0.0 ? -(x) : (x))

double compare_matrices(int m, int n, double **a, int lda, double **b, int ldb)
{
   int i, j;
   double max_diff = 0.0, diff;

   for (i = 0; i < m; i++)
   {
      for (j = 0; j < n; j++)
         diff = abs(a[i][j] - b[i][j]);
      max_diff = (diff > max_diff ? diff : max_diff);
   }

   return max_diff;
}
