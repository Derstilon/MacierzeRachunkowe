#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void MY_MMult(double **, double **, double **, int, int, long long *, long long*);
void MY_MCopyBlock(double **, double **, int, int, int, int, int, int);
void MY_MSubstractBlock(double **, double **, double **, int, int, int, int, int);
void MY_MSumBlock(double **, double **, double **, int, int, int, int, int);

long long mult_count_inv = 0;
long long add_sub_count_inv = 0;

void inverse_matrix_inner(double **a, double **a_inv, int size) {
    if (size == 2) {
        double inv_det = 1 / (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        a_inv[0][0] = inv_det * a[1][1];
        a_inv[1][0] = -inv_det * a[1][0];
        a_inv[0][1] = -inv_det * a[0][1];
        a_inv[1][1] = inv_det * a[0][0];
        mult_count_inv += 7;
        add_sub_count_inv += 1;
    } else {
        int i;
        long long tmp_mult, tmp_add;
        int size2 = size / 2;
        double **A11_inv = (double**) malloc(size2 * sizeof(double*));
        double **S22 = (double **)malloc(size2 * sizeof(double *));
        double **S22_inv = (double **)malloc(size2 * sizeof(double *));
        double **A_tmp1 = (double **)malloc(size2 * sizeof(double *));
        double **A_tmp2 = (double **)malloc(size2 * sizeof(double *));
        for (i = 0;i<size2;i++){
            A11_inv[i] = (double*) malloc(size2 * sizeof(double));
            S22[i] = (double *)malloc(size2 * sizeof(double));
            S22_inv[i] = (double *)malloc(size2 * sizeof(double));
            A_tmp1[i] = (double *)malloc(size2 * sizeof(double));
            A_tmp2[i] = (double *)malloc(size2 * sizeof(double));
        }
        /// CODE
        // step 1 - calculate inv(A11)
        MY_MCopyBlock(a, A_tmp1, size2, 0, 0, 0, 0, 0);
        inverse_matrix_inner(A_tmp1, A11_inv, size2);
        //step 2 - calculate S22
        MY_MCopyBlock(a, S22, size2, size2, 0, 0, 0, 0);

        MY_MMult(S22, A11_inv, A_tmp2, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;

        MY_MCopyBlock(a, A_tmp1, size2, 0, size2, 0, 0, 0);
        MY_MMult(A_tmp2, A_tmp1, S22, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;

        MY_MSubstractBlock(a, S22, S22, size2, size2, size2, 0, 0);
        add_sub_count_inv += size2 * size2;

        //step 3 - invert S22
        inverse_matrix_inner(S22, S22_inv, size2);

        //step 4 - have fun
        MY_MMult(A11_inv, A_tmp1, A_tmp2, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;
        MY_MMult(A_tmp2, S22_inv, A_tmp1, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;

        MY_MCopyBlock(S22_inv, a_inv, size2, 0, 0, size2, size2, 0); // a_inv22 - done
        MY_MCopyBlock(A_tmp1, a_inv, size2, 0, 0, 0, size2, 1); // a_inv12 - done

        MY_MMult(A_tmp1, a+size2, A_tmp2, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;
        MY_MMult(A_tmp2, A11_inv, A_tmp1, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;

        MY_MCopyBlock(A_tmp1, a_inv, size2, 0, 0, 0, 0, 0);
        MY_MSumBlock(A_tmp1, A11_inv, a_inv, size2, 0, 0, 0, 0);  // a_inv11 - done

        MY_MCopyBlock(a, A_tmp1, size2, size2, 0, 0, 0, 0);
        MY_MMult(S22_inv, A_tmp1, A_tmp2, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;
        MY_MMult(A_tmp2, A11_inv, A_tmp1, size2, 8, &tmp_mult, &tmp_add);
        mult_count_inv += tmp_mult;
        add_sub_count_inv += tmp_add;
        MY_MCopyBlock(A_tmp1, a_inv, size2, 0, 0, size2, 0, 1); // a_inv21 - done

        /// CODE
        for (i = 0; i < size2; i++) {
          free(A11_inv[i]);
          free(S22[i]);
          free(S22_inv[i]);
          free(A_tmp1[i]);
          free(A_tmp2[i]);
        }
        free(A11_inv);
        free(S22);
        free(S22_inv);
        free(A_tmp1);
        free(A_tmp2);
    }
}

void inverse_matrix(double **a, double **a_inv, int size, long long *num_mult, long long *num_add) {
    mult_count_inv = 0;
    add_sub_count_inv = 0;
    inverse_matrix_inner(a, a_inv, size);
    *num_mult = mult_count_inv;
    *num_add = add_sub_count_inv;
}