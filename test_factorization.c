#include <stdio.h>
// #include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"

void MY_MMult(double **, double **, double **, int, int, long long *, long long *);
void copy_matrix(int, int, double **, int, double **, int);
void random_matrix(int, int, double **, int);
double compare_matrices(int, int, double **, int, double **, int);
void print_matrix(int, int, double **, int);
double dclock();
void inverse_matrix(double **, double **, int, long long *, long long *);
void lu_factorization(double **, double **, double **, int, long long *, long long *);

double fabs(double x)
{
   if (x > 0)
      return x;
   else
      return -x;
}

int main()
{
   int len, i, j;
   double **a, **a_prim, **l, **u;
   double start_time, end_time, best_time;
   long long mult_count, add_count, tmp_mult_count, tmp_add_count;
   int status = 1;

   for (len = 2; len <= pow(2, LLAST); len *= 2)
   {
      a = (double **)malloc(len * sizeof(double *));
      a_prim = (double **)malloc(len * sizeof(double *));
      l = (double **)malloc(len * sizeof(double *));
      u = (double **)malloc(len * sizeof(double *));
      for (i = 0; i < len; i++)
      {
         a[i] = (double *)malloc(len * sizeof(double));
         l[i] = (double *)malloc(len * sizeof(double));
         u[i] = (double *)malloc(len * sizeof(double));
         a_prim[i] = (double *)malloc(len * sizeof(double));
      }

      random_matrix(len, len, a, 0);

      for (i = 0; i < NREPEATS; i++)
      {
         start_time = dclock();
         lu_factorization(a, l, u, len, &mult_count, &add_count);
         end_time = dclock();
         if (i == 0)
         {
            best_time = end_time - start_time;
         }
         else
         {
            best_time = (end_time - start_time < best_time ? end_time - start_time : best_time);
         }
      }

      /// Error check
      // print_matrix(l, l, a, l);
      // print_matrix(l, l, a_inv, l);

      MY_MMult(l, u, a_prim, len, 8, &tmp_mult_count, &tmp_add_count);
      // print_matrix(l, l, b, l);
      status = 1;

      status = 1;
      if (compare_matrices(len, len, a, 0, a_prim, 0) > EPS)
         status = 0;

      printf("%d %d %le %lld %lld\n", len, status, best_time, mult_count / NREPEATS, add_count / NREPEATS);
      fflush(stdout);

      for (i = 0; i < len; i++)
      {
         free(a[i]);
         free(l[i]);
         free(u[i]);
         free(a_prim[i]);
      }
      free(a);
      free(l);
      free(u);
      free(a_prim);
   }
   return 0;
}