#include <stdio.h>
// #include <malloc.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"

void REF_MMult(int, int, int, double **, int, double **, int, double **, int );
void MY_MMult(double **, double **, double **, int, int, long long *, long long*);
void copy_matrix(int, int, double **, int, double **, int );
void random_matrix(int, int, double **, int);
double compare_matrices( int, int, double **, int, double **, int );
void print_matrix(int m, int n, double **a, int lda);
double dclock();

int main()
{
  int 
    l, 
    m, n, k,
    lda, ldb, ldc, 
    rep;

  double
    dtime, dtime_best,        
    gflops, 
    diff;

  double 
    **a, **b, **c, **cref, **cold;

  long long mult_count, add_count;

  printf("MY_MMult = [\n");

  for ( l=2; l<=pow(2, LLAST); l *=2 ){
    m = l;
    n = l;
    k = l;

    lda = m;
    ldb = n;
    ldc = m;

    /* Allocate space for the matrices */
    /* Note: I create an extra column in A to make sure that
       prefetching beyond the matrix does not cause a segfault */
    a = (double **)malloc(lda * sizeof(double *));
    for (int i = 0; i < k; i++) {
      a[i] = (double *)malloc(k * sizeof(double));
    }
    //a = ( double ** ) malloc( lda * sizeof( double *) );
    b = (double **)malloc(ldb * sizeof(double *));
    for (int i = 0; i < m; i++) {
      b[i] = (double *)malloc(m * sizeof(double));
    }
    //b = ( double ** ) malloc( ldb * sizeof( double *) );
    c = (double **)malloc(ldc * sizeof(double *));
    for (int i = 0; i < m; i++) {
      c[i] = (double *)malloc(m * sizeof(double));
    }
    //c = ( double ** ) malloc( ldc * sizeof( double *) );
    //cold = ( double ** ) malloc( ldc * sizeof( double *) );
    cold = (double **)malloc(ldc * sizeof(double *));
    for (int i = 0; i < m; i++) {
      cold[i] = (double *)malloc(m * sizeof(double));
    }
    //cref = ( double ** ) malloc( ldc * sizeof( double *) );
    cref = (double **)malloc(ldc * sizeof(double *));
    for (int i = 0; i < m; i++) {
      cref[i] = (double *)malloc(m * sizeof(double));
    }

    /* Generate random matrices A, B, Cold */
    random_matrix( m, k, a, lda );
    random_matrix( k, n, b, ldb );
    random_matrix( m, n, cold, ldc );

    copy_matrix( m, n, cold, ldc, cref, ldc );

    /* Run the reference implementation so the answers can be compared */

    REF_MMult( m, n, k, a, lda, b, ldb, cref, ldc );

    /* Time the "optimized" implementation */
    for ( rep=0; rep<NREPEATS; rep++ ){
      copy_matrix( m, n, cold, ldc, c, ldc );

      /* Time your implementation */
      dtime = dclock();

      MY_MMult( a, b, c, l, THRESHOLD, &mult_count, &add_count);
      
      dtime = dclock() - dtime;

      if ( rep==0 )
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }
    //print_matrix(m, n, a, lda);
    //print_matrix(m, n, b, lda);
    //print_matrix(m, n, c, lda);
    //print_matrix(m, n, cref, lda);

    diff = compare_matrices( m, n, c, ldc, cref, ldc );

    printf("%d %le %le %ld %ld \n", l, dtime_best, diff, mult_count, add_count);
    fflush( stdout );

    for (int i = 0; i < lda; i++) {
      free(a[i]);
    }
    free(a);
    for (int i = 0; i < ldb; i++) {
      free(b[i]);
    }
    free(b);
    for (int i = 0; i < ldc; i++) {
      free(c[i]);
    }
    free(c);
    for (int i = 0; i < ldc; i++) {
      free(cold[i]);
    }
    free(cold);
    for (int i = 0; i < ldc; i++) {
      free(cref[i]);
    }
    free(cref);
  }

  printf( "];\n" );

  exit( 0 );
}

