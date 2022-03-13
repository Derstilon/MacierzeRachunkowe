#include <stdlib.h>

int mult_count = 0;
int add_sub_count = 0;

void MY_MSumBlock(double **a, double **b, double **c, int blockSize, int ia, int ja, int ib, int jb)
{
  int k, l;
  for (k = 0; k < blockSize; k++)
  {
    for (l = 0; l < blockSize; l++)
    {
      c[k][l] = a[ia + k][ja + l] + b[ib + k][jb + l];
    }
  }
  add_sub_count += blockSize * blockSize;
}

void MY_MSubstractBlock(double **a, double **b, double **c, int blockSize, int ia, int ja, int ib, int jb)
{
  int k, l;
  for (k = 0; k < blockSize; k++)
  {
    for (l = 0; l < blockSize; l++)
    {
      c[k][l] = a[ia + k][ja + l] - b[ib + k][jb + l];
    }
  }
  add_sub_count += blockSize * blockSize;
}

void MY_MCopyBlock(double **a, double **b, int blockSize, int ia, int ja,
                   int ib, int jb) {
  int k, l;
  for (k = 0; k < blockSize; k++)
  {
    for (l = 0; l < blockSize; l++)
    {
      b[ib + k][jb + l] = a[ia + k][ja + l];
    }
  }
}

void MY_MMultBlockBinet(double **a, double **b, double **c, int blockSize, int ia, int ja, int ib, int jb)
{
  double ***Pa, ***Pb, ***Pc;
  int k, l;
  int blockSize2 = blockSize / 2;
  /* multiply a block of size blockSize x blockSize of a and b and store the result in c */
  if (blockSize == 2)
  {
    c[0][0] = a[ia][ja] * b[ib][jb] + a[ia][ja+1] * b[ib+1][jb];
    c[0][1] = a[ia][ja] * b[ib][jb + 1] + a[ia][ja+1] * b[ib+1][jb+1];
    c[1][0] = a[ia + 1][ja] * b[ib][jb] + a[ia+1][ja+1] * b[ib+1][jb];
    c[1][1] = a[ia + 1][ja] * b[ib][jb+1] + a[ia + 1][ja + 1] * b[ib+1][jb+1];
    mult_count += 8;
    add_sub_count += 4;
  }
  else
  {
    // allocate memory for the matrices
    Pa = (double ***)malloc(4 * sizeof(double **));
    Pb = (double ***)malloc(4 * sizeof(double **));
    Pc = (double ***)malloc(4 * sizeof(double **));
    for (k = 0; k < 4; k++)
    {
      Pa[k] = (double **)malloc(4 * sizeof(double *));
      Pb[k] = (double **)malloc(4 * sizeof(double *));
      Pc[k] = (double **)malloc(4 * sizeof(double *));
      for (l = 0; l < 4; l++)
      {
        Pa[k][l] = (double *)malloc(blockSize2 * sizeof(double));
        Pb[k][l] = (double *)malloc(blockSize2 * sizeof(double));
        Pc[k][l] = (double *)malloc(blockSize2 * sizeof(double));
      }
    }
    /* recursively call the function */

    //(A11*B11) + (A12*B21)
    MY_MMultBlockBinet(a, b, Pa[0], blockSize2, ia, ja, ib, jb);
    MY_MMultBlockBinet(a, b, Pb[0], blockSize2, ia, ja + blockSize2, ib + blockSize2, jb);

    //(A11*B21) + (A12*B22)
    MY_MMultBlockBinet(a, b, Pa[1], blockSize2, ia, ja, ib + blockSize2, jb);
    MY_MMultBlockBinet(a, b, Pb[1], blockSize2, ia, ja + blockSize2, ib + blockSize2, jb + blockSize2);

    //(A21*B11) + (A22*B21)
    MY_MMultBlockBinet(a, b, Pa[2], blockSize2, ia + blockSize2, ja, ib, jb);
    MY_MMultBlockBinet(a, b, Pb[2], blockSize2, ia + blockSize2, ja + blockSize2, ib + blockSize2, jb);

    //(A21*B12) + (A22*B22)
    MY_MMultBlockBinet(a, b, Pa[3], blockSize2, ia + blockSize2, ja, ib, jb + blockSize2);
    MY_MMultBlockBinet(a, b, Pb[3], blockSize2, ia + blockSize2, ja + blockSize2, ib + blockSize2, jb + blockSize2);

    // Calculate block values
    for (k = 0; k < 4; k++)
    {
      MY_MSumBlock(Pa[k], Pb[k], Pc[k], blockSize2, 0, 0, 0, 0);
    }

    /*
     * express matrix c in terms of Pk
     * c=
     * | Pc[0] | Pc[1] |
     * | Pc[2] | Pc[3] |
     */
    MY_MCopyBlock(Pc[0], c, blockSize2, 0, 0, 0, 0);
    MY_MCopyBlock(Pc[1], c, blockSize2, 0, 0, 0, blockSize2);
    MY_MCopyBlock(Pc[2], c, blockSize2, 0, 0, blockSize2, 0);
    MY_MCopyBlock(Pc[3], c, blockSize2, 0, 0, blockSize2, blockSize2);
    for (k = 0; k < 4; k++)
    {
      for (l = 0; l < 4; l++)
      {
        free(Pa[k][l]);
        free(Pb[k][l]);
        free(Pc[k][l]);
      }
      free(Pa[k]);
      free(Pb[k]);
      free(Pc[k]);
    }
    free(Pa);
    free(Pb);
    free(Pc);
  }
}

void MY_MMultBlockStrassen(double **a, double **b, double **c, int blockSize, int i, int j, int threshold)
{
  /* multiply a block of size blockSize x blockSize of a and b and store the result in c */
  double ***Pa, ***Pb, ***Pc;
  int k, l;
  int blockSize2 = blockSize / 2;
  if (blockSize <= threshold)
  {
    MY_MMultBlockBinet(a, b, c, blockSize, 0, 0, 0, 0);
  }
  else
  {
    // allocate memory for the matrices
    Pa = (double ***)malloc(7 * sizeof(double **));
    Pb = (double ***)malloc(7 * sizeof(double **));
    Pc = (double ***)malloc(7 * sizeof(double **));
    for (k = 0; k < 7; k++)
    {
      Pa[k] = (double **)malloc(blockSize2 * sizeof(double *));
      Pb[k] = (double **)malloc(blockSize2 * sizeof(double *));
      Pc[k] = (double **)malloc(blockSize2 * sizeof(double *));
      for (l = 0; l < blockSize2; l++)
      {
        Pa[k][l] = (double *)malloc(blockSize2 * sizeof(double));
        Pb[k][l] = (double *)malloc(blockSize2 * sizeof(double));
        Pc[k][l] = (double *)malloc(blockSize2 * sizeof(double));
      }
    }
    // copy Strassen matrices

    //(A11+A22) * (B11+B22)
    MY_MSumBlock(a, a, Pa[0], blockSize2, i, j, i + blockSize2, j + blockSize2);
    MY_MSumBlock(b, b, Pb[0], blockSize2, i, j, i + blockSize2, j + blockSize2);

    //(A21+A22) * B11
    MY_MSumBlock(a, a, Pa[1], blockSize2, i + blockSize2, j, i + blockSize2, j + blockSize2);
    MY_MCopyBlock(b, Pb[1], blockSize2, i, j, 0, 0);

    // A11 * (B12-B22)
    MY_MCopyBlock(a, Pa[2], blockSize2, i, j + blockSize2, 0, 0);
    MY_MSubstractBlock(b, b, Pb[2], blockSize2, i, j + blockSize2, i + blockSize2, j + blockSize2);

    // A22 * (B21-B11)
    MY_MCopyBlock(a, Pa[3], blockSize2, i + blockSize2, j + blockSize2, 0, 0);
    MY_MSubstractBlock(b, b, Pb[3], blockSize2, i + blockSize2, j, i, j);

    //(A11+A12) * B22
    MY_MSumBlock(a, a, Pa[4], blockSize2, i, j, i, j + blockSize2);
    MY_MCopyBlock(b, Pb[4], blockSize2, i + blockSize2, j + blockSize2, 0, 0);

    //(A21-A11) * (B11+B12)
    MY_MSubstractBlock(a, a, Pa[5], blockSize2, i + blockSize2, j, i, j);
    MY_MSumBlock(b, b, Pb[5], blockSize2, i, j, i, j + blockSize2);

    //(A12-A22) * (B21+B22)
    MY_MSubstractBlock(a, a, Pa[6], blockSize2, i, j + blockSize2, i + blockSize2, j + blockSize2);
    MY_MSumBlock(b, b, Pb[6], blockSize2, i + blockSize2, j, i + blockSize2, j + blockSize2);

    /* recursively call the function */
    for (k = 0; k < 7; k++)
      MY_MMultBlockStrassen(Pa[k], Pb[k], Pc[k], blockSize2, i, j, threshold);

    /*
     * express matrix c in terms of Pk
     * c=
     * | Pc[0]+Pc[3]-Pc[4]+Pc[6] |      Pc[2]+Pc[4]        |
     * |      Pc[1]+Pc[4]        | Pc[0]-Pc[1]+Pc[2]+Pc[5] |
     */
    for (k = 0; k < blockSize2; k++)
    {
      for (l = 0; l < blockSize2; l++)
      {
        c[i + k][j + l] = Pc[0][k][l] + Pc[3][k][l] - Pc[4][k][l] + Pc[6][k][l];
        c[i + k][j + l + blockSize2] = Pc[2][k][l] + Pc[4][k][l];
        c[i + k + blockSize2][j + l] = Pc[1][k][l] + Pc[3][k][l];
        c[i + k + blockSize2][j + l + blockSize2] = Pc[0][k][l] - Pc[1][k][l] + Pc[2][k][l] + Pc[5][k][l];
      }
    }
    add_sub_count += 8 * blockSize2 * blockSize2;
    for (k = 0; k < 7; k++)
    {
      for (l = 0; l < blockSize2; l++)
      {
        free(Pa[k][l]);
        free(Pb[k][l]);
        free(Pc[k][l]);
      }
      free(Pa[k]);
      free(Pb[k]);
      free(Pc[k]);
    }
    free(Pa);
    free(Pb);
    free(Pc);
  }
}

void MY_MMult(double **a, double **b, double **c, int blockSize, int threshold, long long* num_mult, long long* num_add)
{
  MY_MMultBlockStrassen(a, b, c, blockSize, 0, 0, threshold);
  *num_mult = mult_count;
  *num_add = add_sub_count;
}
