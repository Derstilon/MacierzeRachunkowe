#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MERGE_OPERATION_COUNT   \
   mult_count_lu += tmp_mult;   \
   add_sub_count_lu += tmp_add; \
   tmp_mult = 0;                \
   tmp_add = 0

#define INV_BLOCK(m_tmp, m_inv) \
   inverse_matrix(m_tmp, m_inv, size2, &tmp_mult, &tmp_add)

#define INVERSE_MATRIX(matrix)            \
   INV_BLOCK(matrix##_tmp, matrix##_inv); \
   MERGE_OPERATION_COUNT

#define MULTIPLY_MATRIX(a, b, ab)                     \
   MY_MMult(a, b, ab, size2, 8, &tmp_mult, &tmp_add); \
   MERGE_OPERATION_COUNT

#define COPY_BLOCK(source, destination, i, j, k, l, negate) \
   MY_MCopyBlock(source, destination, size2, i *size2, j *size2, k *size2, l *size2, negate);

#define COPY_MATRIX(source, destination, i, j) \
   COPY_BLOCK(source, destination, i, j, 0, 0, 0)

#define PASTE_MATRIX(source, destination, i, j) \
   COPY_BLOCK(source, destination, 0, 0, i, j, 0)

#define SUBSTRACT_MATRIX(a, i, j, b, c)                      \
   MY_MSubstractBlock(a, b, c, i *size2, j *size2, 0, 0, 0); \
   add_sub_count_lu += size2 * size2

#define LU_FACTOR(matrix) \
   lu_factorization(matrix, l_tmp, u_tmp, size2, &tmp_mult, &tmp_add)

void MY_MMult(double **, double **, double **, int, int, long long *, long long *);
void MY_MCopyBlock(double **, double **, int, int, int, int, int, int);
void MY_MSubstractBlock(double **, double **, double **, int, int, int, int, int);
void MY_MSumBlock(double **, double **, double **, int, int, int, int, int);
void inverse_matrix(double **, double **, int, long long *, long long *);

long long mult_count_lu = 0;
long long add_sub_count_lu = 0;

void lu_factorization_inner(double **a, double **l, double **u, int size)
{
   if (size == 2)
   {
      l[0][0] = 1;
      l[0][1] = 0;
      l[1][0] = a[1][0] / a[0][0];
      l[1][1] = 1;
      u[0][0] = a[0][0];
      u[0][1] = a[0][1];
      u[1][0] = 0;
      u[1][1] = a[1][1] - l[1][0] * a[0][1];
      mult_count_lu += 2;
      add_sub_count_lu += 1;
   }
   else
   {
      int i;
      long long tmp_mult = 0, tmp_add = 0;
      int size2 = size / 2;

      double **u_inv = (double **)malloc(size * sizeof(double *));
      double **l_inv = (double **)malloc(size * sizeof(double *));

      double **a_tmp = (double **)malloc(size * sizeof(double *));
      double **l_tmp = (double **)malloc(size * sizeof(double *));
      double **u_tmp = (double **)malloc(size * sizeof(double *));
      double **u_inv_tmp = (double **)malloc(size * sizeof(double *));
      double **l_inv_tmp = (double **)malloc(size * sizeof(double *));
      double **s_tmp1 = (double **)malloc(size * sizeof(double *));
      double **s_tmp2 = (double **)malloc(size * sizeof(double *));

      for (i = 0; i < size; i++)
      {
         u_inv[i] = (double *)malloc(size * sizeof(double));
         l_inv[i] = (double *)malloc(size * sizeof(double));
         a_tmp[i] = (double *)malloc(size * sizeof(double));
         l_tmp[i] = (double *)malloc(size * sizeof(double));
         u_tmp[i] = (double *)malloc(size * sizeof(double));
         u_inv_tmp[i] = (double *)malloc(size * sizeof(double));
         l_inv_tmp[i] = (double *)malloc(size * sizeof(double));
         s_tmp1[i] = (double *)malloc(size * sizeof(double));
         s_tmp2[i] = (double *)malloc(size * sizeof(double));
      }
      /*************** L11 | U11 ***************/
      LU_FACTOR(a);
      PASTE_MATRIX(l, l_tmp, 0, 0);
      PASTE_MATRIX(u, u_tmp, 0, 0);

      /****************** L21 ******************/
      INVERSE_MATRIX(u);

      COPY_MATRIX(a, a_tmp, 1, 0);
      COPY_MATRIX(u_inv, u_inv_tmp, 0, 0);
      MULTIPLY_MATRIX(a_tmp, u_inv, l_tmp);
      PASTE_MATRIX(l_tmp, l, 1, 0);

      /****************** U12 ******************/
      COPY_MATRIX(l, l_tmp, 0, 0);
      INVERSE_MATRIX(l);

      COPY_MATRIX(a, a_tmp, 0, 1);
      COPY_MATRIX(l_inv, l_inv_tmp, 0, 0);
      MULTIPLY_MATRIX(l_inv, a_tmp, u_tmp);
      PASTE_MATRIX(u_tmp, u, 0, 1);

      /******************* S *******************/
      COPY_MATRIX(a, a_tmp, 1, 0);
      COPY_MATRIX(u_inv, u_inv_tmp, 0, 0);
      MULTIPLY_MATRIX(a_tmp, u_inv, s_tmp1);
      COPY_MATRIX(l_inv, l_inv_tmp, 0, 0);
      MULTIPLY_MATRIX(s_tmp1, l_inv, s_tmp2);
      COPY_MATRIX(a, a_tmp, 0, 1);
      MULTIPLY_MATRIX(s_tmp2, a_tmp, s_tmp1);
      SUBSTRACT_MATRIX(a, 1, 1, s_tmp1, s_tmp2);

      /*************** L22 | U22 ***************/
      LU_FACTOR(s_tmp2);
      PASTE_MATRIX(u_tmp, u, 1, 1);
      PASTE_MATRIX(l_tmp, l, 1, 1);

      for (i = 0; i < size; i++)
      {
         free(u_inv[i]);
         free(l_inv[i]);
         free(a_tmp[i]);
         free(l_tmp[i]);
         free(u_tmp[i]);
         free(u_inv_tmp[i]);
         free(l_inv_tmp[i]);
         free(s_tmp1[i]);
         free(s_tmp2[i]);
      }
      free(u_inv);
      free(l_inv);
      free(a_tmp);
      free(l_tmp);
      free(u_tmp);
      free(u_inv_tmp);
      free(l_inv_tmp);
      free(s_tmp1);
      free(s_tmp2);
   }
}

void lu_factorization(double **a, double **l, double **u, int size, long long *num_mult, long long *num_add)
{
   lu_factorization_inner(a, l, u, size);
   *num_mult = mult_count_lu;
   *num_add = add_sub_count_lu;
}