void MY_MMultBlock(double **a, double **b, double **c, int blockSize, int i, int j)
{
  /* multiply a block of size blockSize x blockSize of a and b and store the result in c */
  if (blockSize == 1)
  {
    c[i][j] = a[i][j] * b[i][j];
  }
  else
  {
    /* recursively call the function */
    MY_MMultBlock(a, b, c, blockSize / 2, i, j);
    MY_MMultBlock(a, b, c, blockSize / 2, i, j + blockSize / 2);
    MY_MMultBlock(a, b, c, blockSize / 2, i + blockSize / 2, j);
    MY_MMultBlock(a, b, c, blockSize / 2, i + blockSize / 2, j + blockSize / 2);
  }
}

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
}

void MyMCopyBlock(double **a, double **b, int blockSize, int ia, int ja, int ib, int jb)
{
  int k, l;
  for (k = 0; k < blockSize; k++)
  {
    for (l = 0; l < blockSize; l++)
    {
      b[ib + k][jb + l] = a[ia + k][ja + l];
    }
  }
}

void MY_MMultStrassen(double **a, double **b, double **c, int blockSize, int i, int j, int threshold)
{
  /* multiply a block of size blockSize x blockSize of a and b and store the result in c */
  double ***Pa, ***Pb, ***Pc;
  int k, l;
  if (blockSize <= threshold)
  {
    MY_MMultBlock(a, b, c, blockSize, i, j);
  }
  else
  {
    // allocate memory for the matrices
    Pa = (double ***)malloc(7 * sizeof(double **));
    Pb = (double ***)malloc(7 * sizeof(double **));
    Pc = (double ***)malloc(7 * sizeof(double **));
    for (k = 0; k < 7; k++)
    {
      Pa[k] = (double **)malloc(blockSize * sizeof(double *));
      Pb[k] = (double **)malloc(blockSize * sizeof(double *));
      Pc[k] = (double **)malloc(blockSize * sizeof(double *));
      for (l = 0; l < blockSize; l++)
      {
        Pa[k][l] = (double *)malloc(blockSize * sizeof(double));
        Pb[k][l] = (double *)malloc(blockSize * sizeof(double));
        Pc[k][l] = (double *)malloc(blockSize * sizeof(double));
      }
    }
    // copy Strassen matrices

    //(A11+A22) * (B11+B22)
    MY_MSumBlock(a, a, Pa[0], blockSize, i, j, i + blockSize / 2, j + blockSize / 2);
    MY_MSumBlock(b, b, Pb[0], blockSize, i, j, i + blockSize / 2, j + blockSize / 2);

    //(A21+A22) * B11
    MY_MSumBlock(a, a, Pa[1], blockSize, i + blockSize / 2, j, i + blockSize / 2, j + blockSize / 2);
    MY_MCopyBlock(b, Pb[1], blockSize, i, j, 0, 0);

    // A11 * (B12-B22)
    MY_MCopyBlock(a, Pa[2], blockSize, i, j + blockSize / 2, 0, 0);
    MY_MSubstractBlock(b, b, Pb[2], blockSize, i, j + blockSize / 2, i + blockSize / 2, j + blockSize / 2);

    // A22 * (B21-B11)
    MY_MCopyBlock(a, Pa[3], blockSize, i + blockSize / 2, j + blockSize / 2, 0, 0);
    MY_MSubstractBlock(b, b, Pb[3], blockSize, i + blockSize / 2, j, i, j);

    //(A11+A12) * B22
    MY_MSumBlock(a, a, Pa[4], blockSize, i, j, i, j + blockSize / 2);
    MY_MCopyBlock(b, Pb[4], blockSize, i + blockSize / 2, j + blockSize / 2, 0, 0);

    //(A21-A11) * (B11+B12)
    MY_MSubstractBlock(a, a, Pa[5], blockSize, i + blockSize / 2, j, i, j);
    MY_MSumBlock(b, b, Pb[5], blockSize, i, j, i, j + blockSize / 2);

    //(A12-A22) * (B21+B22)
    MY_MSubstractBlock(a, a, Pa[6], blockSize, i, j + blockSize / 2, i + blockSize / 2, j + blockSize / 2);
    MY_MSumBlock(b, b, Pb[6], blockSize, i + blockSize / 2, j, i + blockSize / 2, j + blockSize / 2);

    /* recursively call the function */
    for (k = 0; k < 7; k++)
      MY_MMultStrassen(Pa[k], Pb[k], Pc[k], blockSize / 2, i, j, threshold);

    // express matrix c in terms of Pk
    for (k = 0; k < blockSize; k++)
    {
      for (l = 0; l < blockSize; l++)
      {
        c[i + k][j + l] = Pc[0][k][l] + Pc[3][k][l] - Pc[4][k][l] + Pc[6][k][l];
        c[i + k + blockSize / 2][j + l] = Pc[1][k][l] + Pc[3][k][l];
        c[i + k][j + l + blockSize / 2] = Pc[2][k][l] + Pc[4][k][l];
        c[i + k + blockSize / 2][j + l + blockSize / 2] = Pc[0][k][l] - Pc[1][k][l] + Pc[2][k][l] + Pc[5][k][l];
      }
    }
  }
}
