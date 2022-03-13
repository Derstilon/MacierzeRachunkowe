/* Create macros so that the matrices are stored in column-major order */

/* Routine for computing C = A * B + C */

void REF_MMult( int m, int n, int k, double **a, int lda, 
                                    double **b, int ldb,
                                    double **c, int ldc )
{
  int i, j, p;

  for ( i=0; i<lda; i++ ){
    for ( j=0; j<ldb; j++ ){
      for ( p=0; p<ldc; p++ ){
	      c[i][j] +=  a[i][p] * b[p][j];
      }
    }
  }
}


  
